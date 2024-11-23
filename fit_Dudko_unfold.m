function [theta,theta_std,resnorm,resid] = fit_Dudko_unfold(pd_obs,edges,Tmean,Fdot,theta0,par)
% Fitting parameters for the DHS or CHS unfolding model
%   DHS: Dudko, Hummer & Szabo (2006), CHS: Cossio, Hummer & Szabo (2016)
% Input:
%   pd_obs   column array of experiment probability densities
%   edges    bin edges for pd_obs.  pd_obs(i) is the mean probabiltydensity
%                  in the interval )bin) from edges(i) to edges(i+1)
%       pd_obs and edges, may be caclulated by probdens.m
%   Tmean        Temperature (°C)
%   Fdot     Mean value of df/dt before unfolding
%   theta0   Initial guess for model parameters [DG;dx;log10(k0)]
%   par      parameter struct with fields Fdot, nu, model
%     nu:  1/2: quadratic potentiial, 2/3 cubic potential
%     model alternatives:
%       'DHS' (Dodko, Hummer & Szabo, 2006)  
%       'CHS' (Cossio, Hummer & Szabo, 2016) 
% Output:
%   theta     Fitted parameters
%   theta_std Standard deviations for the parameter estimates 
%   resnorm   Norm of difference between input and model pd_obs

% Version 4.1,October 2023

  % Optimization parameters:
  opt = optimoptions('lsqcurvefit');
  opt.MaxFunctionEvaluations = 2000;
  opt.Display = 'off';

  F  = (edges(1:end-1)+edges(2:end))'/2;  % Bin midpoints
  
  % theta0: initial guess for theta = [DG;dx;lgk0]
  DG_0 = theta0(1);
  dx_0 = theta0(2);
  log10k0_0 = theta0(3);

  a_0 = par.nu*dx_0/DG_0;  % Use a as parameter.  Setting limits on a
                           % lets us avoid complex numbers in kfun.
  % Alternative parameter vector for use in fminsearch:
  thetacalc0 = [DG_0;a_0;log10k0_0];
  % Symmetric error bars in log10(k) avoids nonsensical error bars
  
  % Lower and upper parameter bounds
  lb = [0;0;-25];  %If theta(3) = log10(k)
  ub = [500;0.9999/F(end);1];  % Make sure that 1-a*F > 0 in the probalility expression
  
  % Function for calculating model probabilities
  probfun = @(thetacalc,F)Dudko_unfold_probability(thetacalc,F,Tmean,Fdot,par);
  
  % Fit model parameters to data:
  [thetacalc,resnorm,resid,exitflag,~,~,J] = ...
    lsqcurvefit(probfun,thetacalc0,F,pd_obs,lb,ub,opt);
  if exitflag < 1
    % warning('lsqcurvefit problems. Exitflag: %d',exitflag)
    theta = NaN;
    theta_std = NaN;
    resnorm = NaN;
    return
  end
  DG = thetacalc(1);
  a = thetacalc(2);
  lgk0 = thetacalc(3);
  dx  = a*DG/par.nu;
  theta = [DG;dx;lgk0];
  
  % ***  Confidence interval for identified perameters ***
  % J is the Jacbian of f(thetacalc) = taufun(thetacalc,...).  
  % JG is the Jacbian of fG(theta) = f(thetacalc(theta))
  da_dDG = -par.nu*dx/DG^2;
  da_dx = par.nu/DG;
  JG = [J(:,1)+J(:,2)*da_dDG,J(:,2)*da_dx,J(:,3)];
  ci = nlparci(theta,resid,'jacobian',JG);  % 95% confidence interval
  theta_std = (ci(:,2)-theta)/fzero(@(x)normcdf(x)-0.975,-1); % Standard deviation
end