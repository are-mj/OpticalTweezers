function [Trip,Tzip,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_experiment(file,plotting,par)
% Generalised version of analyse_experiment that can identify multiple rips 
% per trace and also late rips )i,e, rips in relaxing trace). 
% To identfy multiple rips: set par.maxrips > 1
% To identify late rips: set par.laterips = 1
% To identify late rips in Trip use logical array laterips = Trip.Fdot<0 

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
  else % Called from app
    stack = dbstack;
    folder = datafolder;
    if ~isscalar(stack) & strcmp(stack(2).file,'RipAnalysis.mlapp')
      load RipAnalysis_settings appsettings
      folder = appsettings.Datafolder;
    end
    filename = fullfile(folder,file);
    if isfile(filename)
      shortname = shorten_filename(filename);
    else
      error("File %s not found",file);
    end
  end

  [t,f,x,T] = read_experiment_file(filename);
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
  [t,f,x,T] = remove_time_loops(t,f,x,T);

  
  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  if isempty(peakpos) || isempty(valleypos)
    fprintf("No peaks or valleys in Filename: %s\n",shortname)
    return
  end
  % Decimate time series if number of points per trace is too high:
  [~,~,t,f,x,T] = decim(peakpos,valleypos,par,t,f,x,T);

  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
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
    rngpull = valleypos(end)+1:peakpos(end);  % Bug fix 20250907
    x(rngpull) = x(rngpull)-min(x(rngpull));   
  end  

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
    nzp = numel(r.force);
    if nzp < 1
      r.work = [];
    end
    for zpno = 1:nzp
      r.work (zpno,1) = Crooks_work(r.force(zpno),r.deltax(zpno),...
        r.temperature(zpno),par);
    end
    r = trim_trace(r,par);
    if ~isempty(r.force)  
      relax = [relax;r];
      Tzip = [Tzip;create_table(r)];
    end
  end

  for cycleno = 1:numel(valleypos)-1
    rng = valleypos(cycleno):valleypos(cycleno+1);
    if sum(peakpos>valleypos(cycleno) & peakpos<valleypos(cycleno+1)) ~= 1
      continue  % keep only cycles with exactly one peak
    end
    [~,pkpos] = max(x(rng));
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
    nrp = length(p.force);
    p.cycleno = repmat(cycleno,[nrp,1]);
    p.topforce = repmat(p.f(end),[nrp,1]); 
    p.pullingspeed = abs(median(diff(p.x)./diff(p.t)))*ones(nrp,1); 
    if par.laterips == 0
      bad = isempty(p.force) || p.force < 0;
      % if all(bad)  & Bug, because all([]) is true!
      if ~isempty(bad) & all(bad) & par.laterips == 0% Corected 2025-09-05
        continue  % Skip this cycle
      end
    end
    
    fn = string(fieldnames(p));
    for i = 6:length(fn)
      p.(fn(i))(bad,:) = [];
    end
    nrp = length(p.force);
    if nrp < 1
      p.work = [];  % Bug fix 2025-09-05
    end      
    for rpno = 1:nrp
      p.work(rpno,1) = Crooks_work(p.force(rpno),p.deltax(rpno), ...
        p.temperature(rpno),par);
    end   

    % Relaxation trace struct
    % rlxrng = rng(pkpos+1:end);  % Test 20260119 to include early laterips
    rlxrng = rng(pkpos:end);
    r.file = file;
    r.t = t(rlxrng);
    r.f = f(rlxrng);
    r.x = x(rlxrng);
    r.T = T(rlxrng);    
    r = rip_finder(r,par);
    nzp = length(r.ripx);
    if nzp > 0 & r.force > 0
      r.pullingspeed = abs(median(diff(r.x)./diff(r.t)))*ones(nzp,1);
      for zpno = 1:nzp
        r.work(zpno,1) = Crooks_work(r.force(zpno),r.deltax(zpno), ...
          r.temperature(zpno),par);
      end
      r.topforce = r.f(1)*ones(nzp,1);
      r.cycleno = cycleno*ones(nzp,1);
      r = trim_trace(r,par);
      if ~isempty(r.force)
        relax = [relax;r];
        Tzip = [Tzip;create_table(r)];
      end
    end

    if par.laterips  
      p = laterip_trace(r,p,par);
      p.cycleno = cycleno;
    end
    if isempty(p.force)
      continue
    end
    p = trim_trace(p,par);
    if ~isempty(p.force)
      pull = [pull;p];
      Trip = [Trip;create_table(p)];
    end    

  end
  % Handle pulling trace after last cycle:
  if peakpos(end) > valleypos(end)
    pullrng = valleypos(end):peakpos(end);
    if length(pullrng) > 40  % Eliminate unrealisputically short ranges
      p.file = file;
      p.t = t(pullrng);
      p.f = f(pullrng);
      p.x = x(pullrng);
      p.T = T(pullrng);   
      p.pullingspeed = 0;
      p.cycleno = length(cycleno)+1;
      p.topforce = p.f(end);
      p = rip_finder(p,par);
      nrp = length(p.force);
      p.pullingspeed = abs(median(diff(p.x)./diff(p.t)))*ones(nrp,1);
      p.topforce = p.topforce*ones(length(p.force),1);
      p.cycleno = p.cycleno*ones(length(p.force),1);  
      nrp = length(p.ripx);
      if exist('r',"var")
        nzp = length(r.ripx);
        for zpno = 1:nzp
          r.work(zpno,1) = Crooks_work(r.force(zpno),r.deltax(zpno), ...
            r.temperature(zpno),par);
        end
        for rpno = 1:nrp
          p.work(rpno,1) = Crooks_work(p.force(rpno),p.deltax(rpno),p.temperature(rpno),par);
        end
      else
        p.work = NaN;
      end
      if par.laterips  
        p = laterip_trace(r,p,par);
      end      
      p = trim_trace(p,par);
      if ~isempty(p.force)
        pull = [pull;p];
        Trip = [Trip;create_table(p)];
      end
    end
  end
  if height(Trip)+height(Tzip) < 1
    fprintf("No rips or zips found. Filename: %s\n",shortname);
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
   
    legend(legtext(textelements));
  end
end

function st = trim_trace(st,par)
% Remove structs with unlikely variable values from trace struct 
  bad = false;
  laterip = st.topforce - st.f(1) > 5 & st.topforce - st.f(end) > 5;
  % st comprises both pulling and relaxing trace -> late rip found
  if st.fdot > 0 || laterip % Pulling trace
    deltaxlim = par.deltaxlimits_rips;
  else  % Relaxing trace
    deltaxlim = par.deltaxlimits_zips;
    bad = bad | st.force > st.topforce*par.maxzipfactor;
  end
  bad = bad | st.deltax < deltaxlim(1) | st.deltax > deltaxlim(2);
  if sum(bad)> 0
    fn = string(fieldnames(st));
    % The first five fields of st are not duplicated, so we start at 6
    for i = 6:length(fn)
      st.(fn(i))(bad,:) = [];
    end
  end
end
  