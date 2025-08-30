function [theta,resnorm] = fit_Bell_unfold(pd_obs,edges,Tmean,Fdot,theta0)
% Fitting Bell_unfold_probability to observed probability densities
% Input:
%   pd_obs   column array of experiment probability densities
%   edges    bin edges for pd_obs.  pd_obs(i) is the mean probabiltydensity
%                  in the interval )bin) from edges(i) to edges(i+1)
%       pd_obs and edges, may be caclulated by probdens.m
%   Tmean    Temperature (°C)
%   Fdot     mean value of dF/dt before unfolding 
%   theta0   Initial guess for model parameters [dx;log10(k0)]
% Output:
%   theta     Fitted parameters
%   resnorm   Norm of difference between input and model pd_obs

% Version 2.0 2025-06-25 Removed confidence interval by nlparci
%                        as this often gave incorrect results 

  % Optimization parameters:
  opt = optimoptions('lsqcurvefit');
  opt.MaxFunctionEvaluations = 2000;
  opt.Display = 'off';

  F  = (edges(1:end-1)+edges(2:end))'/2;  % Bin midpoints
  lb = [0;-15];
  ub = [500;10000];

  probfun = @(theta,F) Bell_unfold_probability(theta,F,Tmean,Fdot);
  [theta,resnorm,resid,exitflag,~,~,J] = ...
    lsqcurvefit(probfun,theta0,F,pd_obs,lb,ub,opt);
  if exitflag < 1
    warning('lsqcurvefit problems. Exitflag: %d',exitflag)
  end 
end
