function [Trip,Tzip,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_experiment(file,plotting,par)
% Find rips and zips data in protein stretching experiemnt records
% Input:
%   file:     data file from optical tweezer instrument  
%   plotting: 1 for force vs time plot.  Default: no plot
%   par:      parameter struct.  Default: params
%  Output:
%  Trip, Tzip:  Matlab tables of properties of rips and zips
%  pull, relax: struct arrays with detailed data for pulling and relax
%               traces with an identified rip or zip
%  t,f,x,T:     arrays of time, force, trap position and temperature for
%               complete file
%  peakpos:     record index for force peaks 
%  valleypos:   record indicee fr force valleys

  % Make sure all output variables are defined:
  Trip = [];
  Tzip = [];
  pull = [];
  relax = [];
  t = [];
  f = [];
  x = [];
  T = [];

  if nargin < 3
    par = params;
  end
  if nargin < 2
    plotting = 0;
  end

  file = string(strrep(file,'\','/'));  % Use Unix separator
  data_folder = string(strrep(datafolder,'\','/'));  
  n = strlength(data_folder);
  if isfile(file)  % file contains full path
    filename = file;
    if extractBetween(file,1,n)== data_folder
      % First part of file == datafolder
      shortname = extractAfter(file,n+1);
    else
      shortname = file; % Use full  path
    end
  elseif isfile(fullfile(datafolder,file))
    filename = file;
    shortname = file;
  else
    error("File %s not found",file);
  end

  [t,f,x,T] = read_experiment_file(filename);
  if numel(t)<10
    waitfor(msgbox(sprintf('%s contains no useful data',filename)));
    error('Bad eperiment file')
  end
  
  % Eliminate obviously faulty records
  bad = isnan(x) | isnan(f) | isnan(t);
  f(bad) = []; x(bad)=[]; t(bad) = [];
  [t,f,x,T] = remove_time_loops(t,f,x,T);

  if plotting
    figure;
    plot(t,f);
    hold on;
    xlabel('Time (s)');
    ylabel('Force (pN)')
    title(shortname,'interpreter','none');
  end

  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  if isempty(peakpos) || isempty(valleypos)
    return
  end
  
  % analyse all traces for rips/zips
  npeaks = numel(peakpos);
  nvalleys = numel(valleypos);

  valleyfirst = peakpos(1)>valleypos(1);
  peakfirst = peakpos(1)<valleypos(1);
 
  % Remove drift in x by forcing x=0 at valleys
  % NOTE: This inevitably results in discontinuos x at peaks
  x(1:valleypos(1)) = x(1:valleypos(1))-min(x(1:valleypos(1)));
  for i = 1:numel(valleypos)-1
    rngpull = valleypos(i)+1:peakpos(i+peakfirst);
    x(rngpull) = x(rngpull)-min(x(rngpull));
    rngrlx = peakpos(i+peakfirst)+1:valleypos(i+1);
    x(rngrlx) = x(rngrlx) - min(x(rngrlx));
  end
  if peakpos(end)>valleypos(end)
    rngpull = valleypos(end):peakpos(end);
    x(rngpull) = x(rngpull)-min(x(rngpull));   
  end

  pull = [];  % Array of pulling trace structs
  relax = []; % Array of relaxing trace structs
  Trip = [];    % Pulling results table
  Tzip = [];    % Relaxing results table
  if plotting
    plot(t(peakpos),f(peakpos),'.k',t(valleypos),f(valleypos),'.k');
  end

  peakfirst = valleypos(1)>peakpos(1);
  for i = 1:npeaks-peakfirst
    rng = valleypos(i):peakpos(i+peakfirst);
    if numel(rng)<50  % not a real trace
      continue
    end
    % Pulling trace struct p
    p.t = t(rng);
    p.f = f(rng);
    p.x = x(rng);
    p.T = T(rng);
    p.file = shortname;
    p = singlerip_finder(p,par);
    p.pullingspeed = median(diff(p.x)./diff(p.t));
    % Rips with very low deltax are not likely to be real:
    if ~isempty(p.force) && p.deltax > 5 && p.deltax < 30
      pull = [pull;p];
      Trip = [Trip;create_table(p)];
    end
  end
  for i = 1:nvalleys-1
    rng = peakpos(i):valleypos(i+1-peakfirst);
    if numel(rng)<50
      continue  % not a real trace
    end
    % Relaxing trace struct r
    r.t = t(rng);
    r.f = f(rng);
    r.x = x(rng);
    r.T = T(rng);
    r.file = shortname;
    r = singlerip_finder(r,par);
    r.pullingspeed = abs(median(diff(r.x)./diff(r.t)));
    if ~isempty(r.force) && r.force < f(peakpos(i))*par.maxzipfactor ...
        && r.deltax < -2
      relax = [relax;r];
      Tzip = [Tzip;create_table(r)];
    end 
  end
  if plotting
    if ~isempty(Trip)
      plot(Trip.Time,Trip.Force,'*r')
    end
    if ~isempty(Tzip)
      plot(Tzip.Time,Tzip.Force,'ok',MarkerFaceColor='w');
    end
  end
end