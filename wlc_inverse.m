function x = wlc_inverse(f,P,T,L0,simple)
% Worm-like chain extension as function of the force
% Input:
%  F: Force (pN)
%  P: Persistence length (nm)
%  T: Temperature (K)
%  L0: Contour length (Molecule length at maximun extension) (nm)
%  simple: Do not use the improved fit to WLC (default)
% Output: 
%  x: Extension of molecule (nm) 

% Version 1.1: Makes sure that 0 < x < L0
  
  if nargin < 5
    simple = 1;
  end
  x = zeros(size(f));
  for i = 1:length(f)
    fun = @(x) wlc(x,P,T,L0,simple)-f(i);
    x(i) = fzero(fun,[0,L0*0.9999]);
  end
end