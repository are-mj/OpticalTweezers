function p = Bell_unfold_probability(theta,F,Tmean,Fdotmean)
% Unfold probability density using the Bell model
% Input:
%  theta     : [dx;log10(k0)];
%  Tmean     : Temperature (°C)  (scalar)
%  F         : Force vector (pN)
%  Fdotmean  : Observed rate of change for force

  dx = theta(1);
  lgk0 = theta(2);

  kB = 0.01380649; % Boltzmann constant, unit: 1e-21 J = 1zJ/K
  beta = 1./(kB*(Tmean+273.15));
  k0 = 10^lgk0;

  % Unfolding rate:
  kU = k0*exp(beta*F*dx);
  p = kU/Fdotmean.*exp(k0/(beta*Fdotmean*dx)*(1-exp(beta*F*dx)));
end