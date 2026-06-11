function s = rip_finder(data,range,par)
% Identifies the most likely single  rip/zip in an optical tweezers 
% pulling or relaxing trace. Returns the the input struct with 
% additinal fields specifying details for the event
% Inputs;:
%   data: [t,f,x,T] array of all time series.  (T is oprional)
%   range: two element array of start and end row in data
%   par:  parameter struct (default: params)
%
% Output
%  s: output trace struct with fields:
  %  range
  %  force 
  %  deltax
  %  time
  %  ripx
  %  fdot
  %  fstep
  %  rip_index  
  %  pfx_b     
  %  pfx_a  
  %  dt   
  %  temperature
  %  noise
  %  fitrange
  %  topforce
  %  pullingspeed

% 20260307:  Revised version to handle multiple rips/zips in the same trace

% Author: Are Mjaavatten

  s = empty_trace;
  % s = [];
  if nargin < 3
    par = params;
  end
  if range(2)-range(1) < par.minpointspertrace
    return  % Too short trace
  end
  f = data(range(1):range(2),2);
  % if min(f)<0
  %   return  % Unrealistic data
  % end
  sgn = sign(f(end)-f(1));  % +1 for pull, -1 for relax

  % rng = valid_trace_part(f,sgn,par)+range(1)-1;  % Absolute range
  rng = range(1):range(end);
  if isempty(rng)
    return  % More than two flat parts. Probably nonsense data
  end
  range = [rng(1),rng(end)];
  t = data(rng,1);
  f = data(rng,2);
  x = data(rng,3);
  T_OK = size(data,2) == 4;   
  if T_OK
    T = data(rng,4);
  end  
  n_points = numel(f);
  if n_points<par.minpointspertrace
    return
  end
  smoothwindow = n_points*0.02;
  stdwindow = round(n_points*0.1);

  noise = movstd(f-smoothdata(f,1,'movmean',smoothwindow),stdwindow);
  mean_noise = mean(noise);
  slope = movingslope(f,par.supportlength);
  dslope = detrend(slope);  
  sslope = sort(-dslope,'descend');
  min_peak = max(sslope(10),par.minpeak_slope);  
  
  % Unfolding events give peaks in dslope:
  warning('off','signal:findpeaks:largeMinPeakHeight');
  [~,rip_index] = findpeaks(-sgn*dslope,"MinPeakHeight",min_peak, ...
    "MinPeakDistance",par.supportlength);  
  % rip_index(j) is the index of a point very near the steepest slope 
  % during rip no j.
  warning('on','signal:findpeaks:largeMinPeakHeight');
  if isempty(rip_index)  % No unfoldings found
    return
  end  
  topforce = (sgn>0)*f(end) + (sgn<0)*f(1);  % Start of pull or end of relax trace
  range = [rng(1),rng(end)];
  [~,~,fdot,fstep,weight] = rip_finder_fit(data,range,rip_index,par);
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
  

  % Remove lower quality rips in close rip clusters :
  epsilon = numel(f)/20;
  quality = fstep.*weight;
  quality = quality(valid);
  rip_index = merge_rip_clusters(rip_index,quality,epsilon);
  
  [pfx_b,pfx_a,fdot,fstep,weight,noise,fitrange] = ...
    rip_finder_fit(data,range,rip_index,par);
  
  quality = fstep.*weight;
  [~,ix] = sort(quality,'descend');
  maxrips = min(length(rip_index),par.maxrips);    

  for j = 1:maxrips
    k = ix(j);  % process in descending quality order
    % Search for largest change in f near rip_index(k)
    searchstart = max(rip_index(k)-2*par.supportlength,1);
    searchend = min(rip_index(k)+par.supportlength,n_points);
    [~,pos] = max(diff(-sgn*f(searchstart:searchend)));
    rippos = pos+searchstart-1;
  
    if sgn*fstep(k) >= par.min_fstep  & sgn*quality(k) > 0.35
      s(j,1) = empty_trace;
      s(j,1).range = range;
      s(j,1).ripx = x(rippos);
      s(j,1).pfx_b = pfx_b(k,:);
      s(j,1).pfx_a = pfx_a(k,:);
      % The force is found by linear interpolation at s.ripx
      s(j,1).force = polyval(s(j,1).pfx_b,s(j,1).ripx);
      if sgn > 0  % Treat case where rip force > f(end)
        if s(j,1).force > f(end)
            % Let s.pfx_a be s.pfx_b shifted by fstep
            s(j,1).pfx_a = [s(j,1).pfx_b(1),s(j,1).pfx_b(2)-fstep(k)];
          elseif s(j,1).force < f(1)
            s(j,1) = empty_trace;
            return
        end   
      end
      xend = (s(j,1).force-s(j,1).pfx_a(2))/s(j,1).pfx_a(1);
      s(j,1).deltax = xend - s(j,1).ripx;
      s(j,1).time = t(rippos);
      s(j,1).fdot = fdot(k);
      s(j,1).fstep = fstep(k); 
      s(j,1).rip_index = rip_index(k);
      s(j,1).dt = mean(diff(t));
      if T_OK
        s(j,1).temperature = data(rippos,4);
      else
        s(j,1).temperature = NaN;
      end
      s(j,1).noise = noise(k);
      s(j,1).topforce = topforce;
      s(j,1).fitrange = fitrange(k,:);
      s(j,1).pullingspeed = diff(x(fitrange(k,:)))/diff(t(fitrange(k,:)));
      s(j,1).work = Crooks_work(s(j,1).force,s(j,1).deltax,s(j,1).temperature,par);
    
    end
  end
end


