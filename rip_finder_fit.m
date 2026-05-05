function [pfx_b,pfx_a,fdot,fstep,weight,noise,fitr] ...
  = rip_finder_fit(data,range,rip_index,par)
% Fit linear polynomials to f(x) before and after rips or zips
% Returns data for potentital rips.
% Input: 
%   data:     [t,f,x,T] array of all time series.  (T is oprional)
%   range:    two element array of start and end row in data
%   rip_index: array of rows in daata for potential rips
%   par:  parameter struct (default: params)
% Output: 
%   pfx_b  (n_rips by 2) Polynomial fit to f(x) before rip/zip
%   pfx_a  (n_rips by 2) Polynomial fit to f(x) after rip/zip
%   fdot   (n_rips by 1) Slope of f(t) before rip/zip
%   fstep  (n_rips by 1) value of pfx_b - pfx_a at rip/zip
%   weight (n_rips by 1) weights for ranking rips/zips
%   noise  (n_rips by 1) signal noise before rip/zip
%   fitr   (n_rips by 1)

  if nargin < 4
    par = params;
  end
  t = data(range(1):range(2),1);
  f = data(range(1):range(2),2);
  x = data(range(1):range(2),3);
  n_points = numel(f);      % Number of experiment records
  minpoints = round(n_points*par.minfitfraction);
  maxpoints = round(n_points*par.maxfitfraction);
  % rip_index = [rip_index;n_points-par.ripsteps];  % add virtual rip at end
  n_rips = numel(rip_index);  % NUmber of potential rips or zips
  sgn =  sign(f(end)-f(1));  
  % allocate space
  pfx_b = zeros(n_rips,2);  % Linear polynomial for f(x) before rip
  pfx_a = zeros(n_rips,2);  % Linear polynomial for f(x) after rip
  pft_b = zeros(n_rips,2);  % Linear polynomial for f(t) before rip
  pft_a = zeros(n_rips,2);  % Linear polynomial for f(t) after rip
  noise = zeros(n_rips,1);
  weight = zeros(n_rips,1); 
  fdot  = zeros(n_rips,1);  % Loading rate (df/dt before rip)
  fstep = zeros(n_rips,1);  % rip force change  
  for i = 1:n_rips
    if rip_index(i)<minpoints
      % The potential rip is too close to the start of the trace
      % rip_index(i) = 0;
      continue
    elseif n_points-rip_index(i) < minpoints
      % The potential rip is too close to the end of the trace
      break
    end
    if f(rip_index(i)) > par.overstretch
      continue
    end
    fitstart = max(1,rip_index(i)-maxpoints);    
    fitend = rip_index(i);
    if i < n_rips  % If statement not needed
      fitend = min(fitend,rip_index(i+1)-2);
    end    
    fitrange_b = fitstart:fitend;
    pfx_b(i,:) = polyfit(x(fitrange_b),f(fitrange_b),1);
    pft_b(i,:) = polyfit(t(fitrange_b),f(fitrange_b),1);
    fdot(i) = pft_b(i,1);

    % Fit linear polynomial pfx_a for f(x) after unfolding:
    fitstart = rip_index(i)+1; 
    fitend = min(fitstart + maxpoints,n_points);
    if sgn> 9 && f(rip_index(i)) > f(end)  % use slope from before rip
      f_after = mean(f(fitstart:fitend));
      pfx_a = [pfx_b(i,1),f_after - pfx_a(i,2)*x(rip_index(i))];
      pft_a = [pft_b(i,1),f_after - pft_a(i,2)*x(rip_index(i))];
    else
      fitrange_a = fitstart:fitend;
      fitr(i,:) = [fitrange_b(1),fitrange_a(end)];
      pfx_a(i,:) = polyfit(x(fitrange_a),f(fitrange_a),1);
      pft_a(i,:) = polyfit(t(fitrange_a),f(fitrange_a),1);
    end
    % Test using f(x) to calculate fstep and weight istead of f(t)
    % fstep(i)   = polyval(pft_b(i,:),t(rip_index(i)))-polyval(pft_a(i,:),t(rip_index(i)));
    fstep(i)   = polyval(pfx_b(i,:),x(rip_index(i)))-polyval(pfx_a(i,:),x(rip_index(i)));
    noise(i) = std(f(fitrange_b)-polyval(pft_b(i,:),t(fitrange_b)));

  % Skip events if fstep < noise*par.noisefactors
    noisefactor = (sgn>0)*par.noisefactor(1) + (sgn<0)*par.noisefactor(2);
    ok = sgn*fstep(i) >= max(par.min_fstep,noisefactor*noise(i)) & ...
      pft_b(i,1).*pft_a(i,1)>0;
    % Handle large step changes in x
    ok = ok & max(abs(diff(x(fitr(i,1):fitr(i,2)-1))))<15;  
    % big changes in slope before and after a rip often yields too high fstep
    % Create a weight to counteract this:
    if ok
      % weight(i) = 1./(abs(log(pft_b(i,1)./pft_a(i,1)))+1);
      weight(i) = 1./(abs(log(pfx_b(i,1)./pfx_a(i,1)))+1);
    end      
  end    
end
