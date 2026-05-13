function [Trip,Tzip,data,peakpos,valleypos] ...
  = analyse_experiment(file,plotting,par)
% Find unfolding and refolding evente in optical tweezers experiment file
% Input:  
%   file     - experiment file name
%   plotting - 0: no plot (default), 1: plot force vs. time. Mark evebts.
%   par      - parameter struct.  Dafault: par = params.
% par contains settings and tuning parameters, e.g
%   par.maxrips:  Maximum rips identified per trace (also maximumm zips)
%   par.laterips: Record also rip occuring in the rlaxation trace
%     To identify late rips in TRIP: laterips = Trip.Fdot<0.
% Output:
%   Trip, Tzip: Matlab results tables with one row per rip or zip.  
%     Columns for Rip(or zip) Force, tTme, Deltax and more
%   data: n by 4 array of time series data: data = [t,f,x,T]
%     One row per line in experiment file
%   peakpos, vallepos: row no for force peaks and valleys
%   

  % Make sure all output variables are defined:
  Trip = [];
  Tzip = [];
  pull = [];
  relax = [];
  peakpos = [];
  valleypos = [];

  % Default input:
  if nargin < 3
    par = params;
  end
  if nargin < 2
    plotting = 0;
  end

  % Find experiment file
  file = string(strrep(file,'\','/'));  % Use Unix separator
  data_folder = string(strrep(datafolder,'\','/')); 
  n = strlength(data_folder);
  if isfile(file)  % file contains full path
    filename = file;
    if startsWith(file,data_folder)
      % First part of file == datafolder
      shortname = extractAfter(file,n+1);
    else
      shortname = file; % Use full  path
    end 
  else % file is shortname
    shortname = file;
    folder = datafolder;
    % Check if function is called from the app
    stack = dbstack;
    if ~isscalar(stack) && strcmp(stack(2).file,'RipAnalysis.mlapp')
      load RipAnalysis_settings appsettings
      folder = appsettings.Datafolder;
    end
    filename = fullfile(folder,file);
    if ~isfile(filename)
      error("File %s not found",file);
    end
  end

  data = read_experiment_file(filename);
  if size(data,1) <= 10
    warning('Too few data rows in file')
    return
  end
  f = data(:,2);
  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  if isempty(peakpos) || isempty(valleypos)
    fprintf("No peaks or valleys in Filename: %s\n",shortname)
    return
  end
  % Decimate time series if number of points per trace is too high:
  [data,factor] = decim(data,peakpos,par);
  t = data(:,1);
  f = data(:,2);
  x = data(:,3);
  if size(data,2)> 3
    T = data(:,4);
  end

  if factor > 1  % Repeat peaksandvalleys if time series were decimated
    [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  end
  peakfirst = peakpos(1)<valleypos(1);
  peaklast = peakpos(end)>valleypos(end); 
 
  % Remove drift in x by forcing x=0 at valleys
  % NOTE: This inevitably results in discontinuos x at peaks
  x(1:valleypos(1)) = x(1:valleypos(1))-min(x(1:valleypos(1)));
  for i = 1:numel(valleypos)-1
    rngpull = valleypos(i)+1:peakpos(i+peakfirst);
    x(rngpull) = x(rngpull)-min(x(rngpull));
    rngrlx = peakpos(i+peakfirst)+1:valleypos(i+1);
    x(rngrlx) = x(rngrlx) - min(x(rngrlx));
  end
  % Handle part after last valleypos
  rng2end = valleypos(end)+1:numel(x); 
  x(rng2end) = x(rng2end)-min(x(rng2end));  
  data(:,3) = x;
  
  % First handle relaxing trace before first relax-pull cycle:
  if peakfirst
    r = rip_finder(data,[peakpos(1),valleypos(1)],par);
    for zipno = 1:length(r)
      newr = trim_trace(data,r(zipno),par);
      if ~isempty(newr)
        newr.filename = shortname;
        newr.cycleno = 0;
        relax= newr;
        Trip = create_table(newr);          
      end  
    end    
  end  

  % Full relax-pull traces
  for cycleno = 1:numel(valleypos)-1
    rng = valleypos(cycleno):valleypos(cycleno+1);
    if sum(peakpos>valleypos(cycleno) & peakpos<valleypos(cycleno+1)) ~= 1
      continue  % keep only cycles with exactly one peak
    end
    [~,pkpos] = max(x(rng));
    if pkpos < par.minpointspertrace || ...
        numel(rng)-pkpos < par.minpointspertrace
      continue  % skip very brief traces
    end

    % Check for rips in pulling trace (late rips).  All rips are sorted on
    % fstep and the top par.maxrips rips are retained. If par.laterips == 0
    % a late rip may thus result in lesser rips in the pulling trace being 
    % skipped.
    pullrange = rng([1,pkpos]);
    p = rip_finder(data,pullrange,par);  % Rips in pulling trace
    % Check for late rips. 
    fullrange = [rng(1),rng(end)];
    p_late = laterip(data,fullrange,par);  % Late rips
	% Sort on force shift:
    p1 = [p;p_late];
    [fstep,ix] = sort([p1.fstep],'descend');
    maxrips = min(par.maxrips,length(fstep));
    for ripno = 1:maxrips        
      newp = trim_trace(data,p1(ix(ripno)),par);
      if ~isempty(newp.force)
        if  ~par.laterips
          if newp.fdot < 0  % late rip: skip
            continue
          end
        end
        newp.filename = shortname;
        newp.cycleno = cycleno;
        pull = [pull;newp];
        Trip = [Trip;create_table(newp)];          
      end  
    end  

    % Relaxation trace struct
    relrange = [rng(pkpos+1),valleypos(cycleno+1)];
    r = rip_finder(data,relrange,par);
    [fstep,ix] = sort(-[r.fstep],'descend');
    maxzips = min(par.maxrips,length(fstep));
    for zipno = 1:maxzips
      newr = trim_trace(data,r(ix(zipno)),par);
      if ~isempty(newr.force)
        if newr.force > newr.topforce*par.maxzipfactor
          % Skip zips at high force
          continue;
        end
        newr.filename = shortname;
        newr.cycleno = cycleno;
        relax = [relax;newr];
        Tzip = [Tzip;create_table(newr)];          
      end  
    end 
  end       

  % Handle pulling trace after last full cycle
  if peaklast
    pullrange = [valleypos(end),peakpos(end)];
    p = rip_finder(data,pullrange,par);
    p = trim_trace(data,p,par);
    nrp = length(p);
    for ripno = nrp:-1:1
      bad = isempty(p(ripno).force) || p(ripno).force < 0;
      if bad
        p(ripno) = [];
      end
      p(ripno).filename = shortname;
      p(ripno).cycleno = cycleno+1;  
      pull = [pull;p(ripno)];
      Trip = [Trip;create_table(p(ripno))]; 
    end   
  end

  if height(Trip)+height(Tzip) < 1
    fprintf("No rips or zips found. Filename: %s\n",shortname);
    return
  end
  if plotting
    figure;
    plot(t,f);
    legtext = ["force",'Rip','Late rip','Zip'];
    textelements = 1;
    hold on;
    if ~isempty(Trip)
      plot(Trip.Time(Trip.Fdot>0),Trip.Force(Trip.Fdot>0),'*r');
      textelements = [textelements,2];
      if any(Trip.Fdot<0)
        plot(Trip.Time(Trip.Fdot<0),Trip.Force(Trip.Fdot<0),'^r');
        textelements = [textelements,3];
      end
    end
    if ~isempty(Tzip)
      plot(Tzip.Time,Tzip.Force,'ok','MarkerFaceColor','w');
      textelements = [textelements,4];
    end
    title(shortname,'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Force (pN)');
    legend(legtext(textelements));
    drawnow;
  end
end
  