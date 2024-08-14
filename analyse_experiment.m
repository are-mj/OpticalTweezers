
function [Tp,Tr,pull,relax,t,f,x] = analyse_experiment(file,plotting)
% Find rips and zip data in experiment file

  if nargin < 2
    plotting = 0;
  end
  % Allow 
  if ischar(file)
    if strcmp(file(1:3),'C:\') | strcmp(file(1:3),'C:/') % Full path
      filename = string(file);
    end
  else
    % datafolder = "Clean_files";  % OK for Karina files.  Generalize!!
    filename = fullfile(datafolder,file);
  end
  [t,f,x,T] = read_experiment_file(filename);
  
  % Eliminate obviously faulty records
  bad = isnan(x) | isnan(f) | isnan(t);
  f(bad) = []; x(bad)=[]; t(bad) = [];
  [t,f,x,T] = remove_time_loops(t,f,x,T);
  
  % Make sure no f value == threshold:
  threshold = mean(f);
  exact_hit = f==threshold;
  f(exact_hit) = f(exact_hit)*(1+5*eps);  

  if plotting
      clf;
      plot(t,f);
      hold on;
  end

  % Create arrayy of sampling time orders of magnitude:
  dtmag0 = round(log10(diff(t)));
  % Smooth out brief periods of aberrant smoing times
  % Find borders for sections with identical smoothed orders of madnitude 
  section_borders = rle_smooth(dtmag0,5);
  Tp = [];
  Tr = [];
  pull = [];
  relax = [];
  for k = 2:numel(section_borders)
    range = section_borders(k-1):section_borders(k);
    tt = t(range);
    ff = f(range);
    xx = x(range);
    TT = T(range);
    
    [Tp_,Tr_,pull_,relax_] = analyse_section(tt,ff,xx,TT,filename,threshold,plotting);
    Tp = [Tp;Tp_];
    Tr = [Tr;Tr_];
    pull_ = [pull;pull_];
    relax = [relax;relax_];
    if plotting
      xlabel('Time (s)');
      ylabel('Force (pN)')
      title(file)
    end
end

function [Tp,Tr,pull,relax] = analyse_section(t,f,x,T,filename,threshold,plotting)

  Tp = [];
  Tr = [];
  pull = [];
  relax = [];
  
  % Find all points where f crosses the thresold value:
  up = find(f(2:end)-threshold>0 & f(1:end-1)-threshold<0);
  down = find(f(2:end)-threshold<0 & f(1:end-1)-threshold>0);
 
  nup = numel(up);
  ndown = numel(down);
  if nup > 0
    upfirst = up(1)<down(1); 
  end
  upmax = min(nup,ndown);
  downmax = min(ndown,nup);
  if upmax+downmax <1
    return
  end

  % Eliminate close crossings due to noise
  if upfirst
    k = find(down(1:downmax)-up(1:upmax) < max(down(1:downmax)-up(1:upmax))/10);
    if ~isempty(k)
      up(k) = [];
      down(k) = [];
    end
  else
    k = find(up(1:upmax)-down(1:downmax) < max(up(1:upmax)-down(1:downmax))/10);
    if ~isempty(k)
      up(k) = [];
      down(k) = [];
    end
  end

  % Redefine nup and ndown, without noise-relat3ed crossings
  nup = numel(up);
  ndown = numel(down);
  upfirst = up(1)<down(1);
  uplast = up(end)>down(end);
  downlast = ~uplast;

  % Find valleys and peaks:
  peakpos = [];
  valleypos = [];
  if upfirst
    % Denote the lowest point between 1 and up(1) as a valley pooint
    [~,m] = min(f(1:up(1)));
    valleypos = m;    
    for i = 1:ndown
      [~,m] = max(f(up(i):down(i+1-upfirst)));
      peakpos = [peakpos;up(i)+m-1]; 
    end
    for i = 1:nup-1
      [~,m] = min(f(down(i):up(i+1)));
      valleypos = [valleypos;down(i)+m-1];
    end
  else
    % denote the highest point beween 1 and down(1) as  a peak point
    [~,m] = max(f(1:down(1)));
    peakpos = m;
    for i = 1:ndown-downlast
      [~,m] = min(f(down(i):up(i)));
      valleypos = [valleypos;down(i)+m-1]; 
    end    
    for i = 1:nup-uplast
      [~,m] = max(f(up(i):down(i+1)));
      peakpos = [peakpos;up(i)+m-1];
    end  
  end
  % Denote last extremum point as peak or valley:
  if uplast
    [~,m]=max(f(up(end):end));
    peakpos = [peakpos;up(end)+m-1];
  elseif downlast
    [~,m]=min(f(down(end):end));
    valleypos = [valleypos;down(end)+m-1];    
  end

  % analyse all traces for rips/zips
  npeaks = numel(peakpos);
  nvalleys = numel(valleypos);

  if (f(peakpos(1))-f(valleypos(1)))*(x(peakpos(1))-x(valleypos(1))) < 0
    x = max(x)-x;   % So f and x have same phase
  end
  
  pull = [];  % Array of pulling trace structs
  relax = []; % Array of relaxing trace structs
  Tp = [];    % Pulling results table
  Tr = [];    % Relaxing results table
  if plotting
    plot(t(peakpos),f(peakpos),'.k',t(valleypos),f(valleypos),'.k');
  end

  % for i = 1:npeaks
  %   text(t(peakpos(i)),f(peakpos(i))+0.2,num2str(i));
  % end


  peakfirst = valleypos(1)>peakpos(1);
  for i = 1:npeaks-peakfirst
    s.t = t(valleypos(i):peakpos(i+peakfirst));
    s.f = f(valleypos(i):peakpos(i+peakfirst));
    s.x = x(valleypos(i):peakpos(i+peakfirst));
    s.T = T(valleypos(i):peakpos(i+peakfirst));
    s.pullingspeed = median(diff(s.x)./diff(s.t));
    s.file = filename;
    s = singlerip_finder(s,par_single);
    % Rips with very low deltax are not likely to be real
    if ~isempty(s.force) && s.deltax > 5
      if plotting
        plot(s.time,s.force,'*r');
      end
      pull = [pull;s];
      Tp = [Tp;create_table(s)];
    end
  end
  for i = 1:nvalleys-1
    r.t = t(peakpos(i):valleypos(i+1-peakfirst));
    r.f = f(peakpos(i):valleypos(i+1-peakfirst));
    r.x = x(peakpos(i):valleypos(i+1-peakfirst));
    r.T = T(peakpos(i):valleypos(i+1-peakfirst));
    r.pullingspeed = abs(median(diff(r.x)./diff(r.t)));
    r.file = filename;
    r = singlerip_finder(r,par_single);
    r.file = filename;
    if ~isempty(r.force) && r.force < f(peakpos(i))*0.33 && r.deltax < -2
      if plotting
        plot(r.time,r.force,'ok',MarkerFaceColor='w');
      end
      relax = [relax;r];
      Tr = [Tr;create_table(r)];
    end 
  end
