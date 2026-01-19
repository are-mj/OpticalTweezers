function p = laterip_trace(s,p,par)
% Looks for rips in relaxation trace and adds adds the largest to the 
% pulling trace struct p
% Input s is the relaxing trace. Only s.t s.f, s.x  and s.T are used.
% Output:
%  p: trace data for full pull_relax cycle
%   If valid rips are found, data for the largest rip are used

  sgn = sign(s.f(end)-s.f(1));  % +1 for pull, -1 for relax
  if sgn > 0
    warning('laterip requires the first input to be a relaxing trace struct')
    return
  end
  f = s.f;
  x = s.x;
  t = s.t;  
  n_points = numel(f);
  if n_points<50
    return
  end

  slope = movingslope(f,max(round(n_points*0.01),2));
  dslope = detrend(slope);

  sslope = sort(dslope,'ascend');
  min_peak = max(sslope(10),par.minpeak_slope);
  % Unfolding events give peaks in -slope:
  warning('off','signal:findpeaks:largeMinPeakHeight');
  % rip_index(i) is the index of a point very near the steepest slope 
  % during rip no i.
  [peaks,rip_index] = findpeaks(-dslope,"MinPeakHeight",min_peak, ...
    "MinPeakDistance",par.supportlength);  
  warning('on','signal:findpeaks:largeMinPeakHeight');
  % trip = t(rip_index);
  if isempty(rip_index)  % No unfoldings found
    return
  end

  [~,order] = sort(peaks,'descend');
  rip_index = rip_index(order(1));  % Seect only the largest
  if rip_index < 3
    return
  end

  % minpoints = round(n_points*par.minfitfraction);
  maxpoints = round(n_points*par.maxfitfraction);
  i = 1;
  fitstart = max(1,rip_index(i)-maxpoints);
  fitend = rip_index(i)-par.ripsteps;  
  fitrange_b = fitstart:fitend;
  pfx_b(i,:) = polyfit(x(fitrange_b),f(fitrange_b),1);
  pft_b(i,:) = polyfit(t(fitrange_b),f(fitrange_b),1);  
  fitend = min(rip_index(i) + maxpoints,n_points);
  fitrange_a = rip_index(i):fitend;    
  pfx_a(i,:) = polyfit(x(fitrange_a),f(fitrange_a),1);
  pft_a(i,:) = polyfit(t(fitrange_a),f(fitrange_a),1);
  fstep(i)   = polyval(pft_b(i,:),t(rip_index(i)))-polyval(pft_a(i,:),t(rip_index(i)));
  noise(i) = std(f(fitrange_b)-polyval(pft_b(i,:),t(fitrange_b)));   
  % Skip events if fstep < noise*10 (higher limit than in rip_finder)
  ok = fstep >= noise*3;    
  if ok
    s.ripx = s.x(rip_index(i));
    s.force = polyval(pfx_b(i,:),s.ripx);
    xend = (s.force-pfx_a(i,2))/pfx_a(i,1);
    s.deltax = xend - s.ripx;
    s.time = s.t(rip_index(i));
    
    s.fdot = pft_b(1);
    s.fstep = fstep(i);
    s.rip_index = rip_index(i);  
    s.pfx_b = pfx_b(i,:);     
    s.pfx_a = pfx_a(i,:);  
    s.dt    = mean(diff(s.t));
    s.temperature = s.T(rip_index(i));
    s.noise = noise(i);
    s.fitrange = [fitstart,fitend];
    s.work = NaN;
    s.cycleno = p.cycleno;
    s.topforce = p.topforce;
    if isempty(s.force) || s.deltax < 5 || s.deltax > 30 ||s.force<10
      return
    else 
      p= mergestructs(p,s);
    end    
  end
end

function p = mergestructs(p,s)
  if isempty(p.fstep) | s.fstep > p.fstep
    rip_index = length(p.t)+s.rip_index;
    s.pullingspeed = mean(diff(s.x)./diff(s.t));
    p.t = [p.t;s.t];
    p.f = [p.f;s.f];
    p.x = [p.x;s.x];
    p.T = [p.T;s.T];
    fn = string(fieldnames(p));
    
    for i = 5:length(fn)
      p.(fn(i)) = s.(fn(i));
    end
    p.rip_index = rip_index;
    p.topforce = max(p.f);
  end
end
    