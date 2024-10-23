function [Tp,Tr,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_experiment(file,par,plotting)
% Find rips and zips data in protein stretching experiemnt records
% Input:
%   file:     data file from optical tweezer instrument
%   par:      parameter struct.  Default: parameter_struct
%   plotting: 1 for force vs time plot.  Default: no plot
%  Output:
%  Tp, Tr:      Matlab tables of properties of rips (Tp) and zips (Tr)
%  pull, relax: struct arrays with detailed data for pulling and relax
%               traces with an identified rip or zip
%  t,f,x,T:     arrays of time, force, trap position and temperature for
%               complete file
%  peaks, valleys: record index for force peaks and valleys, separating
%               individual traces

  if nargin < 3
    plotting = 0;
  end
  if nargin < 2
    par = parameter_struct;
  end
  % Allow file name containing full path
  if isfile(file)
    filename = file;
  else
    filename = fullfile(datafolder,file);
    if ~isfile(filename)
      error(sprintf("File %s is not found",filename));
    end
  end
  Tlist = NaN;
  if isfield(par,'Tlist')
    Tlist = par.Tlist;
  end
  [t,f,x,T] = read_experiment_file(filename,Tlist);
  
  % Eliminate obviously faulty records
  bad = isnan(x) | isnan(f) | isnan(t);
  f(bad) = []; x(bad)=[]; t(bad) = [];
  [t,f,x,T] = remove_time_loops(t,f,x,T);

  if plotting
      clf;
      plot(t,f);
      hold on;
  end
  if plotting
    xlabel('Time (s)');
    ylabel('Force (pN)')
    title(file,'interpreter','none');
  end

  [peakpos,valleypos] = peaksandvalleys(f,par.threshold,par.lim,0);
  
  % analyse all traces for rips/zips
  npeaks = numel(peakpos);
  nvalleys = numel(valleypos);

  % Check if f and x have the same phase
  % if (f(peakpos(1))-f(valleypos(1)))*(x(peakpos(1))-x(valleypos(1))) < 0
  %   x = max(x)-x;   % So f and x have same phase
  % end  
  % NOT ROBUST in a few cases. Use correation coefficient instead:
  X = corrcoef(f,x);
  if X(2,1) < 0   % x and f have opposide phase
    x = max(x)-x;
  end


  % Remove drift in x by forcing x=0 at valleys
  x(1:valleypos(1))=x(1:valleypos(1))-x(valleypos(1));
  for i = 2:numel(valleypos)
    rng = valleypos(i-1):valleypos(i);
    p = polyfit([rng(1),rng(end)],x([rng(1),rng(end)]),1);  
    x(rng(1:end-1)) = x(rng(1:end-1))-polyval(p,rng(1:end-1))';
  end  

  pull = [];  % Array of pulling trace structs
  relax = []; % Array of relaxing trace structs
  Tp = [];    % Pulling results table
  Tr = [];    % Relaxing results table
  if plotting
    plot(t(peakpos),f(peakpos),'.k',t(valleypos),f(valleypos),'.k');
  end

  % For debugging: Show peak numbering in plot
  % for i = 1:npeaks
  %   text(t(peakpos(i)),f(peakpos(i))+0.2,num2str(i));
  % end

  peakfirst = valleypos(1)>peakpos(1);
  for i = 1:npeaks-peakfirst
    rng = valleypos(i):peakpos(i+peakfirst);
    if numel(rng)<50  % not a real trace
      continue
    end
    s.t = t(rng);
    s.f = f(rng);
    s.x = x(rng);
    s.T = T(rng);
    s.file = filename;
    s = singlerip_finder(s,par);
    s.pullingspeed = median(diff(s.x)./diff(s.t));
    % Rips with very low deltax are not likely to be real:
    if ~isempty(s.force) && s.deltax > 5
      if plotting
        plot(s.time,s.force,'*r');
      end
      pull = [pull;s];
      Tp = [Tp;create_table(s)];
    end
  end
  for i = 1:nvalleys-1
    rng = peakpos(i):valleypos(i+1-peakfirst);
    if numel(rng)<50
      continue  % not a real trace
    end
    r.t = t(rng);
    r.f = f(rng);
    r.x = x(rng);
    r.T = T(rng);
    r.file = filename;
    r = singlerip_finder(r,par);
    r.pullingspeed = abs(median(diff(r.x)./diff(r.t)));
    r.file = filename;
    % plot(r.t,r.f,'r')
    if ~isempty(r.force) && r.force < f(peakpos(i))*0.5 && r.deltax < -2
      if plotting
        plot(r.time,r.force,'ok',MarkerFaceColor='w');
      end
      relax = [relax;r];
      Tr = [Tr;create_table(r)];
    end 
  end
end