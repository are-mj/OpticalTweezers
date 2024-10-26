function s = singlerip_finder(s,par)
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
%   slope       % mean df/dx before event
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

  s.force = []; 
  s.deltax = [];
  s.time = [];
  s.ripx = [];
  s.fdot = [];
  s.slope = [];
  s.fstep = [];
  s.rip_index = [];  
  s.pfx_b = [];     
  s.pfx_a = [];  
  s.dt    = [];
  s.temperature = [];
  s.noise = [];


  sgn = sign(s.f(end)-s.f(1));  % +1 for pull, -1 for relax
  
  % Eliminate invalid parts of the trace, such as flat parts at start or
  % end
  rng = valid_trace_part(s.f,sgn);
  s.f = s.f(rng);
  s.x = s.x(rng);
  s.t = s.t(rng);
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
    min_peak = max(sslope(10),0.1);
  else  % Use relax trace parameters. Scale by signal noise.
    % Skip first third of trace (refolding unlikely here)
    slope = movingslope(f,max(round(n_points*0.03),2));
    dslope = detrend(slope);
    noise = movstd(f-smoothdata(f,1,'movmean',smoothwindow),stdwindow);
    mean_noise = mean(noise);
    % sslope = sort(slope,'descend');
    sslope = sort(dslope,'descend');
    min_peak = max(sslope(10),0.003);
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
  [~,~,~,~,~,fstep,weight] = singlerip_finder_fit(s,rip_index,par);
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
  if sum(valid) < 1
    return
  end
  rip_index = rip_index(valid);
  fstep = fstep(valid);
  maxrips = sum(valid);
  weight = weight(valid);
  n_rips = min(par.maxrips,maxrips);
  [~,order] = sort(sgn*fstep.*weight,'descend');
  rip_index = sort(rip_index(order(1:n_rips)));
  
  % Repeat fitting after invalid rips are removed:
  [s.pfx_a,s.pfx_b,~,~,fdot,fstep,weight,noise]= singlerip_finder_fit(s,rip_index,par);
  [~,best] = max(fstep.*weight);
  if weight == 0
    return
  end

  % Search for largest diff(f) near best rip_index
  searchstart = max(rip_index(best)-2*par.supportlength,1);
  searchend = min(rip_index(best)+par.supportlength,n_points);
  [~,pos] = max(diff(-sgn*f(searchstart:searchend)));
  rippos = pos+searchstart-1;

  if sgn*fstep >= par.min_fstep
    s.ripx = x(rippos);
    % The force is found by linear interpolation at s.ripx
    s.force = polyval(s.pfx_b(best,:),s.ripx);
    xend = (s.force-s.pfx_a(best,2))/s.pfx_a(best,1);
    s.deltax = xend - s.ripx;
    s.time = t(rippos);
    s.fdot = fdot(best);
    s.slope = s.pfx_b(best,1);  
    s.fdot = fdot(best);
    s.slope = s.pfx_b(best,1);
    s.fstep = fstep(best); 
    s.rip_index = rip_index(best);
    s.pfx_b = s.pfx_b(best,:);
    s.pfx_a = s.pfx_a(best,:);
    s.dt = mean(diff(s.t));
	  s.temperature = s.T(best);
    s.noise = noise(best);
  end
end
