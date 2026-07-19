function p = laterip(data,fullrange,par) 
  p = empty_trace;
  f = data(:,2);
  x = data(:,3);
  t = data(:,1);
  rng = fullrange(1):fullrange(2);
  if any(f(fullrange)<0)
    return
  end
  % pos: Position within rng
  % row: Row no. in data  
  
  [~,toppos] = max(x(rng));   % Position of peak x in rng
  toprow = toppos + rng(1)-2; % absolute position
  d = round(toppos/3);        % points inspected in relaxing trace  
  laterows = toprow+1:min(toprow+d,rng(end));
  topforce = f(toprow);

% ****  Local indexing in laterows *****
  ff = f(laterows);
  tt = t(laterows);
  xx = x(laterows);
  TT = data(laterows,1);
  oklaterows = valid_trace_part(ff,1,par);
  if isempty(oklaterows)
    return;
  end
  % Local indexing in oklaterows:
  okff = ff(oklaterows);  
  okxx = xx(oklaterows);
  oktt = tt(oklaterows);
  okTT = TT(oklaterows);
  slope = -movingslope(okff,round(max(d/50,2)));
  sslope = sort(slope,'descend');
  min_peak = max(sslope(5),par.minpeak_slope);
  % Unfolding events give peaks in -slope:
  warning('off','signal:findpeaks:largeMinPeakHeight');
  % rip_index is the index of a point very near the steepest slope 
  % during rip no i.
  [peaks,rip_index] = findpeaks(slope,"MinPeakHeight",min_peak, ...
    "MinPeakDistance",par.supportlength);  
  warning('on','signal:findpeaks:largeMinPeakHeight');
  if isempty(rip_index)  % No unfoldings found
    return
  end
  [~,order] = sort(peaks,'descend');
  rip_index = rip_index(order(1));  % Select only the largest
  if rip_index < 3 || oklaterows(end) - rip_index < 3
    return  % Too few points for polyfit
  end
  n_points = length(oklaterows);
  maxpoints = round(n_points*par.maxfitfraction);

  % Indices relative to oklaterows
  fitstart_b = max(1,rip_index-maxpoints);
  fitend_b = rip_index-par.ripsteps;  
  fitrange_b = fitstart_b:fitend_b;
  fitstart_a = rip_index+par.ripsteps;
  fitend_a = min(fitstart_a + maxpoints,n_points);
  fitrange_a = fitstart_a:fitend_a;
  if min(fitend_b-fitstart_b,fitend_a-fitstart_a)<2
    return  % Must have at least three points for fitting
  end

  pft_b = polyfit(oktt(fitrange_b),okff(fitrange_b),1);
  pft_a = polyfit(oktt(fitrange_a),okff(fitrange_a),1);
  pfx_b = polyfit(okxx(fitrange_b),okff(fitrange_b),1);
  pfx_a = polyfit(okxx(fitrange_a),okff(fitrange_a),1);   

  xrip  = okxx(rip_index);
  frip  = polyval(pfx_b,xrip);
  fstep = frip-polyval(pfx_a,xrip);
  noise = std(okff(fitrange_b)-polyval(pft_b,oktt(fitrange_b))); 

  ok = fstep >= max(noise*3,par.min_latefstep);
  if ok
    s = p;
    s.ripx = xrip;
    s.force = frip;
    s.deltax = fstep/pfx_a(1);
    s.time = oktt(rip_index);
    
    s.fdot = pft_b(1);
    s.fstep = fstep;
    s.rip_index = rip_index + oklaterows(1) -1 + toppos;  
    s.pfx_b = pfx_b;     
    s.pfx_a = pfx_a;  
    s.dt    = mean(diff(oktt));
    s.temperature = okTT(rip_index);
    s.noise = noise;
    fitrange = [fitstart_b,fitend_b];
    s.fitrange = fitrange + oklaterows(1) -1 + toppos;
    s.work = Crooks_work(s.force,s.deltax,s.temperature,par);
    s.noise = noise;
    s.topforce = topforce;
    s.pullingspeed = diff(x(fitrange))/diff(t(fitrange));  
      
    if isempty(s.force) || s.deltax < par.deltaxlimits_rips(1) ...
        || s.deltax > par.deltaxlimits_rips(2) ||s.force>par.overstretch
      return
    else 
      p= s;
    end    
  end
end
