function [theta,theta_std,resnorm,resid] = fit_Bell_unfold(pd_obs,edges,Tmean,Fdot,theta0)
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
%   theta_std Standard deviations for the parameter estimate 
%   resnorm   Norm of difference between input and model pd_obs

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
    error('lsqcurvefit problems. Exitflag: %d',exitflag)
  end 
  if nargout > 1
    ci = nlparci(theta,resid,'jacobian',J);  % 95% confidence interval
    theta_std = (ci(:,2)-theta)/fzero(@(x)normcdf(x)-0.975,-1); % Standard deviation
  end
end
