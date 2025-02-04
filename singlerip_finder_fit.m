function [pfx_a,pfx_b,pft_a,pft_b,fdot,fstep,weight,noise,fitr] = singlerip_finder_fit(s,rip_index,par)
% Fit linear polynomials to f(x) before and after rips or zips
% Returns data for potentital rips.
% Input: s: stretching or relaxing trace with experiment results:
%            s.t:  time
%            s.f:  force
%            s.x:  trap position
%        rip_index: rows in s.t etc. for potential rips or zips
%        par:       Tuning parameter struct (optional)
% Output: pfx_b  (n_rips by 2) Polynomial fit to f(x) before rip/zip
%         pfx_a  (n_rips by 2) Polynomial fit to f(x) after rip/zip
%         fdot   (n_rips by 1) Slope of f(t) before rip/zip
%         fstep  (n_rips by 1) value of pfx_b - pfx_a at rip/zip

  f = s.f;
  x = s.x;
  t = s.t;

  n_points = numel(s.f);      % Number of experiment records
  minpoints = round(n_points*par.minfitfraction);
  maxpoints = round(n_points*par.maxfitfraction);
  n_rips = numel(rip_index);  % NUmber of potential rips or zips
  rip_index = [rip_index;n_points-par.ripsteps];  % add virtual rip at end
  sgn =  sign(f(end)-f(1));

  % allocate space
  pfx_b = zeros(n_rips,2);  % Linear polynomial for f(x) before rip
  pfx_a = zeros(n_rips,2);  % Linear polynomial for f(x) after rip
  pft_b = zeros(n_rips,2);  % Linear polynomial for f(t) before rip
  pft_a = zeros(n_rips,2);  % Linear polynomial for f(t) after rip
  noise = zeros(n_rips,1);
  
  fdot  = zeros(n_rips,1);  % Loading rate (df/dt before rip)
  fstep = zeros(n_rips,1);  % rip force change

  % Run through all potential unfoldings to find the force steps
  for i = 1:n_rips
    if rip_index(i)<minpoints
      % The potential rip is too close to the start of the trace
      rip_index(i) = 0;
      continue
    % elseif n_points-rip_index(i)-par.supportlength < minpoints
    elseif n_points-rip_index(i) < minpoints
      % The potential rip is too close to the end of the trace
      break
    end
    if s.f(rip_index(i)) > par.overstretch
      rip_index(i) = 0;
      continue
    end
    % Fit linear polynomial pfx_b for f(x) after unfolding:
    fitstart = max(1,rip_index(i)-maxpoints);
    fitend = rip_index(i)-par.ripsteps;
    if fitend-fitstart < 2  % Too few fitting points
      fstep(i) = 0;
      continue
    end
    if i < n_rips
      fitend = min(fitend,rip_index(i+1)-2);
    end
    fitrange_b = fitstart:fitend;
    pfx_b(i,:) = polyfit(x(fitrange_b),f(fitrange_b),1);
    pft_b(i,:) = polyfit(t(fitrange_b),f(fitrange_b),1);
    fdot(i) = pft_b(1);

    % Fit linear polynomial pfx_a for f(x) after unfolding:
    fitstart = rip_index(i); 
    fitend = min(fitstart + maxpoints,n_points);
    fitrange_a = fitstart:fitend;
    fitr = [fitrange_b(1),fitrange_a(end)];
    pfx_a(i,:) = polyfit(x(fitrange_a),f(fitrange_a),1);
    pft_a(i,:) = polyfit(t(fitrange_a),f(fitrange_a),1);
    fstep(i)   = polyval(pft_b(i,:),t(rip_index(i)))-polyval(pft_a(i,:),t(rip_index(i)));
    noise(i) = std(f(fitrange_b)-polyval(pft_b(i,:),t(fitrange_b)));

    % Useful commands to study the rip/zip:
    % figure; plot(t,f,t(rip_index(i)),f(rip_index(i)),'*r');hold on;
    % plot(t(fitrange_b),polyval(pft_b(i,:),t(fitrange_b)),'k',t(fitrange_a),polyval(pft_a(i,:),t(fitrange_a)),'k')
    % figure; plot(x,f,x(rip_index(i)),f(rip_index(i)),'*r');hold on
	  % plot(x(fitrange_b),polyval(pfx_b(i,:),x(fitrange_b)),'k',x(fitrange_a),polyval(pft_a(i,:),x(fitrange_a)),'k')
  end

  % Skip events if fstep < noise*par.noisefactors
  noisefactor = (sgn>0) * par.noisefactor(1) + (sgn<0) *par.noisefactor(2);
  ok = sgn*fstep >= max(par.min_fstep,noisefactor*noise) & pft_b(:,1).*pft_a(:,1)>0;
  % big changes in slope before and after a rip often yields too high fstep
  % Create a weight to counteract this:
  if any(ok)
    % weight = sqrt(nanmin(slopediff(ok))./slopediff);  
    weight = 1./(abs(log(pft_b(:,1)./pft_a(:,1)))+0.5);
  else
    weight = zeros(n_rips,1);
  end
end
