function W = Crooks_work(force,deltax,T,par)
% Folding or unfolding work as in Fit_Crooks.m
% 20250628: Interpret T < 150 as °C
  P = par.WLC_P;
  L0 = par.WLC_L0;
  if T < 150
    T = T + 273.15;  % deg C to K
  end
  T(isnan(T)) = 290;  % Typical value if T is not given
  deltax = abs(deltax);
  Ws = stretchwork(force,deltax,P,T,L0);
  deltax0 = wlc_inverse(force,P,T,L0);
  W = force.*deltax-deltax/deltax0*Ws;  
end

function W = stretchwork(force,deltax,P,T,L0)
% Calcuate the work done stretching to the unfolding force
%
% Input:
%  force: Unfolding force (nm)
%  P: Persistence length (nm)
%  T: Temperature (K)
%  L0: Contour length (Molecule length at maximun extension) (nm)
% simple: Do not use the improved fit to WLC

  if nargin < 6
    simple = 1;
  end
  x0 = 0;
  W = zeros(size(force));
  for i = 1:numel(force)
    x1 = wlc_inverse(force(i),P,T,L0,simple);
    scale = deltax(i)/x1;
    fun = @(x) wlc(x,P,T,L0,simple);
    W(i) = integral(fun,x0,x1)*scale;
  end
end