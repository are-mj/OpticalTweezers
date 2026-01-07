function s = rip_finder(s,par)
% Identifies the most likely single  rip/zip in an optical tweezers 
% pulling or relaxing trace. Returns the the input struct with 
% additinal fields specifying details for the event
% Inputs;:
%   s: input trace struct with fields:
%   x: extent (trap position)
%   f: force 
%   t: time
% Output
%   s: output trace struct with fields:
%   x: extent (trap position)
%   f: force 
%   t: time
%   force:      % pulling force at rip/zip start
%   deltax:     % rip/zip extent
%   ripx:       % event x value
%   fdot:       % rate of change of f before event
%   fstep       % force shift
%   rip_index   % rows of rips/zips in t,x,f arrays
%   pfx_b       % Linear fit to f(x) before rips/zips
%   pfx_a       % Linear fit to f(x) after rips/zips

% Author: Are Mjaavatten
% Version: 0.1  2023-11-13 Handles both pulling and relaxing traces
%   Better handling of short sequnces (mean(a*x-f)
% Version 0.2 2024_01_21: Moved polynomial fitting to subfunction 
%   rip_finder_fit and repeated fitting after pruning rip candidates
% Version 0.5 2024_08_26 Use valid_trace_part to eliminate irrelevant parts
%    of the trace. Added sampling time step dt to trace struct
% Version 0.6 2025_05_27: Corrected bug: fstep, noise and fitrange were not 
%    reduced to only best value before creating final s.


  s.force = []; 
  s.deltax = [];
  s.time = [];
  s.ripx = [];
  s.fdot = [];
  s.fstep = [];
  s.rip_index = [];  
  s.pfx_b = [];     
  s.pfx_a = [];  
  s.dt    = [];
  s.temperature = [];
  s.noise = [];
  s.fitrange = [];

  sgn = sign(s.f(end)-s.f(1));  % +1 for pull, -1 for relax
  
  % Eliminate invalid parts of the trace, such as flat parts at start or
  % end
  rng = valid_trace_part(s.f,sgn);
  s.f = s.f(rng);
  s.x = s.x(rng);
  s.t = s.t(rng);
  s.T = s.T(rng);
  s = lookforrip(s,sgn,par);
end

function s = lookforrip(s,sgn,par)
  f = s.f;
  x = s.x;
  t = s.t;
  n_points = numel(f);
  if n_points<30
    return
  end

  smoothwindow = n_points*0.02;
  stdwindow = round(n_points*0.1);
  if sgn > 0
    slope = movingslope(f,par.supportlength);
    dslope = detrend(slope);
    mean_noise = std(f-smoothdata(f,1,'movmean',smoothwindow));   
    % sslope = sort(-slope,'descend');
    sslope = sort(-dslope,'descend');
    min_peak = max(sslope(10),par.minpeak_slope);
  else  % Use relax trace parameters. Scale by signal noise.
    % Skip first third of trace (refolding unlikely here)
    % slope = movingslope(f,max(round(n_points*0.03),2));
    slope = movingslope(f,par.supportlength);
    dslope = detrend(slope);
    noise = movstd(f-smoothdata(f,1,'movmean',smoothwindow),stdwindow);
    mean_noise = mean(noise);
    sslope = sort(dslope,'descend');
    min_peak = max(sslope(10),par.minpeak_slope);
  end
  
  % Unfolding events give peaks in -slope:
  warning('off','signal:findpeaks:largeMinPeakHeight');
  % rip_index(i) is the index of a point very near the steepest slope 
  % during rip no i.
  [~,rip_index] = findpeaks(-sgn*dslope,"MinPeakHeight",min_peak, ...
    "MinPeakDistance",par.supportlength);  
  warning('on','signal:findpeaks:largeMinPeakHeight');

  % figure;plot(-sgn*(slope))
  if isempty(rip_index)  % No unfoldings found
    return
  end
  [~,~,~,~,fdot,fstep,weight] = rip_finder_fit(s,rip_index,par);
  if isempty(fstep)
    return
  end
  valid = sgn*fstep.*weight > mean_noise*par.noisefactor((3-sgn)/2);
  if sgn<0
    % zips in first thrid of trace are not realistic
    valid = valid & rip_index > n_points*0.33 &f(rip_index)>min(f)+1;
  else
    % Very early rips are often not realistic
    valid = valid & rip_index > n_points*0.1;
  end
  % eliminate (rare) cases where the slope has the wrong sign
  valid = valid & sgn*fdot > 0;  %
  if sum(valid) < 1
    return
  end
  rip_index = rip_index(valid);
  maxrips = sum(valid);
  n_rips = min(par.maxrips,maxrips);
  
  % Repeat fitting after invalid rips are removed:
  [s.pfx_a,s.pfx_b,~,~,fdot,fstep,weight,noise,s.fitrange] = ...
    rip_finder_fit(s,rip_index,par);

  quality = fstep.*weight;
  epsilon = numel(s.f)/20;
  if n_rips > 1
    % If rips are very close, retain only one per cluster
    okrips = merge_rip_clusters(rip_index,quality,epsilon);
    rip_index = rip_index(okrips);
    s.pfx_b = s.pfx_b(okrips,:);
    s.pfx_a = s.pfx_a(okrips,:);
    s.fitrange = s.fitrange(okrips,:);
  else
    % [~,best] = max(fstep.*weight);
    [~,best] = max(abs(fstep).*weight); % abs necessary for relaxing trace
                           % This bug went undiscovered until Aug. 2025!!
    rip_index = rip_index(best);
    fstep = fstep(best);
    s.pfx_b = s.pfx_b(best,:);
    s.pfx_a = s.pfx_a(best,:);
    fdot = fdot(best);
    noise = noise(best);
    s.fitrange = s.fitrange(best,:);
  end

  for i = 1:length(rip_index)
    % Search for largest diff(-sgn*f) near best rip_index
    searchstart = max(rip_index(i)-2*par.supportlength,1);
    searchend = min(rip_index(i)+par.supportlength,n_points);
    [~,pos] = max(diff(-sgn*f(searchstart:searchend)));
    rippos = pos+searchstart-1;
  
    if sgn*fstep(i) >= par.min_fstep
      % s.file(i,1) = file;
      s.ripx(i,1) = x(rippos);
      % The force is found by linear interpolation at s.ripx
      s.force(i,1) = polyval(s.pfx_b(i,:),s.ripx(i));
      xend = (s.force(i)-s.pfx_a(i,2))/s.pfx_a(i,1);
      s.deltax(i,1) = xend - s.ripx(i);
      s.time(i,1) = t(rippos);
      s.fdot(i,1) = fdot(i);
      s.fstep(i,1) = fstep(i); 
      s.rip_index(i,1) = rip_index(i);
      s.dt(i,1) = mean(diff(s.t));
      if isfield(s,'T')
	      s.temperature(i,1) = s.T(rippos);
      end
      s.noise(i,1) = noise(i);
    end
  end
end


function okrips = merge_rip_clusters(x,quality,epsilon)
    % Perform DBSCAN clustering
    idx = dbscan(x, epsilon, 1);
    
    % Initialize logical array
    okrips = false(size(x));
    
    % Get unique cluster labels
    uniqueClusters = unique(idx(idx > 0));  % Exclude noise points (idx == 0)
    
    % Loop through each cluster to find the highest quality point
    for i = 1:length(uniqueClusters)
        cluster = uniqueClusters(i);
        clusterIdx = find(idx == cluster);
        [~, highestQualityPos] = max(quality(clusterIdx));
        okrips(clusterIdx(highestQualityPos)) = true;
    end
end
