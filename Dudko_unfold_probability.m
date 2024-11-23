function p = Dudko_unfold_probability(thetacalc,F,T,Fdot,par)
% Probability density (dp/dF) for unfolding at given F.
%  thetacalc: parameter vector: [DG;a;lgk0]  (lg: log base 10)
%     NOTE that this is different from theta [DG;dx;log10k]
%  F: Force vector
%  T: Temperature (°C)
%  Fdot:  df/dx (pN/s)
% par: struct of parameters (nu, model)  
%   Models: 'DHS': Dudko, Hummer & Szabo (2006)
%           'CHS': Cossio, Hummer & Szabo (2016)

% Are Mjaavatten
% Version 1.1 October 2023.  Replaces unfold_probability.m

  % Parameters
  DG = thetacalc(1);
  a  = thetacalc(2);
  lgk0 = thetacalc(3);
  k0 = 10^lgk0;
  nu = par.nu;

  kB = 0.01380649; % Boltzmann constant, unit: 1e-21 J = 1zJ/K
  beta = 1./(kB*(T+273.15));

  k = kfun(F,nu,DG,a,k0,beta,par.model);
  if strcmp(par.model,'DHS')
    exponent = exp_ana(F,nu,DG,a,k0,beta,Fdot);
  else
    exponent = zeros(numel(F),1);
    for i = 1:length(F)
      exponent(i) = integral(@(F) -kfun(F,nu,DG,a,k0,beta,par.model),0,F(i))/Fdot;
    end
  end
  p = k/Fdot.*exp(exponent);
end

function k = kfun(F,nu,DG,a,k0,beta,model)
  y = min(a*F,1);
  switch model
    case 'CHS'
      k = k0*(1-y).^(2-1/nu).*exp(beta*DG*(1-(1-y).^(1/nu)));
    case 'DHS'
      k = k0*(1-y).^(1/nu-1).*exp(beta*DG*(1-(1-y).^(1/nu)));
    otherwise
      error('Unknown model name: %s',model);
  end
end

function I = exp_ana(F,nu,DG,a,k0,beta,Fdot)
  % Integral in the exponent expressed analytically
  dx = a*DG/nu;
  I = k0/(dx*Fdot*beta)*(1-exp(beta*DG*(1-(1-nu*F*dx/DG).^(1/nu))));
end
  