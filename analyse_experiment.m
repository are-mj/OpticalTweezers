function [Trip,Tzip,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_experiment(file,plotting,par)
% Generalised version of analyse_experiment that can identify multiple  
% rips per trace and also late rips )i,e, rips in relaxing trace). 
% To identfy multiple rips: set par.maxrips > 1
% To identify late rips: set par.laterips = 1
% To identify late rips in Trip use logical array laterips = Trip.Fdot<0 
% Input:
%   file:     data file from optical tweezer instrument  
%   plotting: 1 for force vs time plot.  Default: no plot
%   par:      parameter struct.  Default: params
%  Output:
%    Trip, Tzip:  Matlab tables of properties of rips and zips
%    pull, relax: st	ruct arrays with detailed data for pulling and relax
%               traces with an identified rip or zip
%    t,f,x,T:     arrays of time, force, trap position and temperature for
%               complete file
%    peakpos:     record index for force peaks 
%    valleypos:   record indicee fr force valleys

  % Make sure all output variables are defined:
  Trip = [];
  Tzip = [];
  pull = [];
  relax = [];

  if nargin < 3
    par = params;
  end
  if nargin < 2
    plotting = 0;
  end

  [t,f,x,T] = read_experiment_file(file);
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
  else 
    if isfile(fullfile(datafolder,file))
      filename = file;
      shortname = file;
    else
      error("File %s not found",file);
    end
  end

  if numel(t)<10
    warning('%s contains no useful data',filename)
    return
  end
  if isempty(T)
    T = NaN(size(t));
  end
  
  % Eliminate obviously faulty records
  bad = isnan(x) | isnan(f) | isnan(t);
  f(bad) = []; x(bad)=[]; t(bad) = [];
  % try
  % [t,f,x,T] = remove_time_loops(t,f,x,T);
  % catch
  % end

  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  % downsample so that record per trace < 1000
  factor = round(mean(diff(peakpos))/2000);
  if factor > 1
    t = downsample(movmean(t,factor),factor);
    f = downsample(movmean(f,factor),factor);
    x = downsample(movmean(x,factor),factor);
    T = downsample(movmean(T,factor),factor);
    [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
    % t = decimate(t,factor);
    % f = decimate(f,factor);
    % x = decimate(x,factor);
    % T = decimate(T,factor);
    % [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);    
  end
  if max(numel(peakpos),numel(valleypos)) < 2
    warning("File must contain more than one cycle");
    return
  end

  % analyse all traces for rips/zips

  % valleyfirst = peakpos(1)>valleypos(1);
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

  % First handle relaxing trace before first cycle:
  if peakfirst
    rlxrng = peakpos(1):valleypos(1);
    r.file = file;
    r.t = t(rlxrng);
    r.f = f(rlxrng);
    r.x = x(rlxrng);
    r.T = T(rlxrng);
    r = rip_finder(r,par);
    nzp = length(r.ripx);
    r.cycleno = 0*ones(nzp,1);
    r.pullingspeed = abs(median(diff(r.x)./diff(r.t)))*ones(nzp,1);
    r.topforce = r.f(1)*ones(nzp,1);
    r.work = Crooks_work(r.force,r.deltax,r.temperature,par);
    r = trim_trace(r,par.deltaxlimits_zips);
    if ~isempty(r.force)
      relax = [relax;r];
      Tzip = [Tzip;create_table(r)];
    end
  end

  for cycleno = 1:numel(valleypos)-1
    rng = valleypos(cycleno):valleypos(cycleno+1);
    [~,pkpos] = max(x(rng));
    % [~,pkpos] = max(f(rng));
    if pkpos < 50 || numel(rng)-pkpos < 50
      continue  % skip very brief traces
    end
    % Rips in pulling trace
    pullrng = rng(1:pkpos);
    p.file = file;
    p.t = t(pullrng);
    p.f = f(pullrng);
    p.x = x(pullrng);
    p.T = T(pullrng); 
    p = rip_finder(p,par);
    nrp = length(p.ripx);
    for ripno = 1:nrp
      if isempty(p.force) | p.force(ripno) < 0
        continue
      end
      p.pullingspeed = repmat(median(diff(p.x)./diff(p.t)),[nrp,1]);
      p.cycleno = repmat(cycleno,[nrp,1]);
      p.topforce = repmat(p.f(end),[nrp,1]); 
      work = Crooks_work(p.force(ripno),p.deltax(ripno),p.temperature(ripno),par);
      p.work = repmat(work,[nrp,1]);
    end

    % Initialise relaxation trace struct
    rlxrng = rng(pkpos+1:end);
    r.file = file;
    r.t = t(rlxrng);
    r.f = f(rlxrng);
    r.x = x(rlxrng);
    r.T = T(rlxrng);    
    r.pullingspeed = 0;  % Placeholder only
    r.cycleno = cycleno;
    r.topforce = r.f(1);
    if par.laterips  
      p = laterip_trace(r,p,par);
    end
    p = trim_trace(p,par.deltaxlimits_rips);
    if ~isempty(p.force)
      pull = [pull;p];
      Trip = [Trip;create_table(p)];
    end
    
    % Zips in relaxing trace
    r = rip_finder(r,par);
    nzp = length(r.ripx);
    if nzp < 1 | r.force < 0
      continue
    end
    r.pullingspeed = abs(median(diff(r.x)./diff(r.t)))*ones(nzp,1);
    r.work = Crooks_work(r.force,r.deltax,r.temperature,par);
    r.topforce = r.f(1)*ones(nzp,1);
    r.cycleno = cycleno*ones(nzp,1);
    if ~isempty(r.force)
      relax = [relax;r];
      Tzip = [Tzip;create_table(r)];
    end

  end
  % Handle pulling trace after last cycle:
  if peakpos(end) > valleypos(end)
    pullrng = valleypos(end):peakpos(end);
    p.file = file;
    p.t = t(pullrng);
    p.f = f(pullrng);
    p.x = x(pullrng);
    p.T = T(pullrng);   
    p.pullingspeed = 0;
    p.cycleno = length(cycleno)+1;
    p = rip_finder(p,par);
    p.topforce = p.f(end);
    p.pullingspeed = median(diff(p.x)./diff(p.t))*ones(length(p.force),1);
    p.topforce = p.topforce*ones(length(p.force),1);
    p.cycleno = (cycleno+1)*ones(length(p.force),1);
    p.work = Crooks_work(p.force,p.deltax,p.temperature,par);
    p = trim_trace(p,par.deltaxlimits_rips);
    if ~isempty(p.force)
      pull = [pull;p];
      Trip = [Trip;create_table(p)];
    end
  end
  if plotting
    figure;
    plot(t,f);
    legtext = ["force",'Rip','Late rip','Zip'];
    textelements = 1;
    hold on;
    if ~isempty(Trip)
      plot(Trip.Time(Trip.Fdot>0),Trip.Force(Trip.Fdot>0),'*r');
      plot(Trip.Time(Trip.Fdot<0),Trip.Force(Trip.Fdot<0),'^r');
      textelements = [textelements,2];
    end
    if ~isempty(Tzip)
      plot(Tzip.Time,Tzip.Force,'ok','MarkerFaceColor','w');
      textelements = [textelements,4];
    end
    title(shortname,'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Force (pN)');
   
    % legend shows only symbols preent in plot:
    legend(legtext(textelements),'Location','best');
  end
end

function st = trim_trace(st,deltaxlim)
% Remove structs from trace struct array st if st.deltax is not within
% limits
  if diff(deltaxlim) > 0
    bad = st.deltax < deltaxlim(1) | st.deltax > deltaxlim(2);
  else
    bad = st.deltax > deltaxlim(1) | st.deltax < deltaxlim(2);
  end
  if sum(bad)> 0
    fn = string(fieldnames(st));
    % The first five fields of st are not duplicated, so we start at 6
    for i = 6:length(fn)
      try
        st.(fn(i))(bad) = [];
      end
    end
  end
end
  