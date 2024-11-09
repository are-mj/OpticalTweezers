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
  
  % Allow file name containing full path
  if isfile(file)
    filename = file;
  else
    filename = fullfile(datafolder,file);
    if ~isfile(filename)
      error("File %p is not found",filename);
    end
  end
  Tlist = NaN;
  if isfield(par,'Tlist')
    Tlist = par.Tlist;
  end
  [t,f,x,T] = read_experiment_file(filename,Tlist);
  if numel(t)<10
    waitfor(msgbox(sprintf('%s contains no useful data',filename)));
    error('Bad eperiment file')
  end
  
  % Eliminate obviously faulty records
  bad = isnan(x) | isnan(f) | isnan(t);
  f(bad) = []; x(bad)=[]; t(bad) = [];
  [t,f,x,T] = remove_time_loops(t,f,x,T);

  if plotting
    clf;
    plot(t,f);
    hold on;
    xlabel('Time (s)');
    ylabel('Force (pN)')
    title(file,'interpreter','none');
  end

  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  if isempty(peakpos) || isempty(valleypos)
    return
  end
  
  % analyse all traces for rips/zips
  npeaks = numel(peakpos);
  nvalleys = numel(valleypos);
 
  % Remove drift in x by forcing x=0 at valleys
  x(1:valleypos(1)) = x(1:valleypos(1))-x(valleypos(1));
  for i = 2:numel(valleypos)
    rng = valleypos(i-1):valleypos(i);
    px = polyfit([rng(1),rng(end)],x([rng(1),rng(end)]),1);  
    x(rng(1:end-1)) = x(rng(1:end-1))-polyval(px,rng(1:end-1))';
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
    p.file = filename;
    p = singlerip_finder(p,par);
    p.pullingspeed = median(diff(p.x)./diff(p.t));
    % Rips with very low deltax are not likely to be real:
    if ~isempty(p.force) && p.deltax > 5
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
    r.file = filename;
    r = singlerip_finder(r,par);
    r.pullingspeed = abs(median(diff(r.x)./diff(r.t)));
    if ~isempty(r.force) && r.force < f(peakpos(i))*0.5 && r.deltax < -2
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