function Tout = run_fit(TRIP,TZIP,par)
% Fits parameters for  Bell_unfold, Bell_refold, Dudko, DGkin and DGCrooks
% and calcuates parameter 95% conficence interval using bootstrapping
% function bootci
% Output: Table of fitted parameters for all sets of themperature and 
% pulling speed 
%
% Supersedes earlier versions that used anaytical confidence intervals  
% from function nlparci. This often gives incorrect results in combiantion
% with bounded parameters search.

% Version 2.0 2025-06-25 Calculate confidence intervals by bootci

  if nargin < 3
    par = params;
  end
  conversion = 0.1439326;  % Energy units kcal/kmol
  thetaBell0Rip = [1;-3];
  thetaBell0Zip = [4;4];
  theta0Dudko = [100;2;-3];
  nboot = 100;   % Bootstrap samples for calculating confdence intervals
  kB = 0.01380649; % Boltzmann's constant (zJ/K/molecule)

  Tout = cell2table(cell(0,15),'VariableNames',{'Text','dx',...
    'cidx','log10k0Zip','cilog10k0'...
    'DGDudko','ciDGDudko','dxDudko','cidxDudko','log10k0Dudko',...
    'cilog10k0Dudko','DGkin','ciDGkin','DGCrooks','ciDGCrooks'});
  [ucases,rcases] = cases(TRIP,TZIP);
  dF = 1;
  k = 1;  % Rowno in Tout
  for i = 1:numel(ucases)  
    if sum(ucases(i).selected) > 10
      Tmean = mean(TRIP.Temperature(ucases(i).selected));
      Fdot = mean(TRIP.Fdot(ucases(i).selected));
      forceZip = TZIP.Force(rcases(i).selected);
      Text = rcases(i).text;
    
      % Refold
      Clusterno = 0;
      DGDudko = 0;
      ciDGDudko = [0,0];
      dxDudko = 0;
      cidxDudko = [0,0];
      log10k0Dudko = 0;
      cilog10k0Dudko = [0,0];
      DGkin = 0;
      ciDGkin = [0,0];
      DGCrooks = 0;
      ciDGCrooks = [0,0];
  
      usel = ucases(i).selected;
      rsel = rcases(i).selected;
      texts = strcat(ucases(i).text,rcases(i).text);
      forceZip = TZIP.Force(rsel);
      thfun = @(f) BellfunZip(f,dF,Tmean,Fdot,thetaBell0Zip);
      thR = thfun(forceZip);
      dx = thR(1);
      log10k0Zip = thR(2);
      ciTHZip = bootci(nboot,thfun,forceZip)';
      cidx = ciTHZip(1,:);
      cilog10k0 = ciTHZip(2,:);
      cilogk0widthZip = kB*(Tmean+273.15)*diff(cilog10k0); % For use in DGkin
      Tout = [Tout;table(Text,dx,cidx,log10k0Zip,cilog10k0,DGDudko,...
        ciDGDudko,dxDudko,cidxDudko,log10k0Dudko,cilog10k0Dudko,DGkin,...
        ciDGkin,DGCrooks,ciDGCrooks)];
      k = k+1;

      % Unfold
      Text = ucases(i).text;
      usel = ucases(i).selected;
      rsel = rcases(i).selected;
      [G(i),Gciu(2*i+[1,2])] = fit_Crooks(TRIP(usel,:),TZIP(rsel,:),par,0);
  
      force = TRIP.Force(usel);
  
      % Bell_unfold
      thfun = @(f) BellfunRip(f,dF,Tmean,Fdot,thetaBell0Rip);
      theta = thfun(force);
      dx = theta(1);
      log10k0Rip = theta(2);
      cith = bootci(nboot,thfun,force)';
      cidx = cith(1,:);
      cilog10k0 = cith(2,:);
  
      % Dudko
      thfun = @(f) Dudkofun(f,dF,Tmean,Fdot,theta0Dudko);
      theta = thfun(force);
      DGDudko = theta(1)*conversion;
      dxDudko = theta(2);
      log10k0Dudko = theta(3);
      try
        cith = bootci(nboot,thfun,force)';
      catch
        cith = bootci(nboot,thfun,force)';
      end
      ciDGDudko = cith(1,:)*conversion;
      cidxDudko = cith(2,:);
      cilog10k0Dudko = cith(3,:);
  
      DGkin = -kB*(Tmean+273.15)*(log10k0Rip-log10k0Zip)*log(10)*conversion;
      cilogk0widthRip = kB*(Tmean+273.15)*diff(cilog10k0);
      ciDGkin = DGkin + sqrt(cilogk0widthRip.^2 + cilogk0widthZip.^2)/2*...
        [-1,1]*log(10)*conversion;
  
      DGCrooks = G(i)*conversion;
      ciDGCrooks = Gciu(2*i+[1,2])*conversion;
  
      Tout = [Tout;table(Text,dx,cidx,log10k0Zip,cilog10k0,DGDudko,...
        ciDGDudko,dxDudko,cidxDudko,log10k0Dudko,cilog10k0Dudko,DGkin,...
        ciDGkin,DGCrooks,ciDGCrooks)];
      k = k+1;
    end
  end
  % Move the refolding lines to the end
  unfoldingrows = find(contains(Tout.Text,"Unfolding"));
  refoldingrows = find(contains(Tout.Text,"Refolding"));
  Tout = Tout([unfoldingrows;refoldingrows],:);
end

function [th,rms] = BellfunRip(fs,dF,Tmean,Fdot,theta0)
% Bell function for pulling trace
  [pd,edges] = probdens(fs,dF);
  [th,rms] = fit_Bell_unfold(pd,edges,Tmean,Fdot,theta0);
end

function [th,rms] = BellfunZip(fs,dF,Tmean,Fdot,theta0)
% Bell function for relaxing trace
  [pd,edges] = probdens(fs,dF);
  [th,rms] = fit_Bell_refold(pd,edges,Tmean,Fdot,theta0);
  if rms > 0.05
    % If the fit is bad, try with another theta0
    theta0 = [9;4];
    [th,rms] = fit_Bell_refold(pd,edges,Tmean,Fdot,theta0);
    if rms > 0.05
      warning('No convergence');
    end
  end
end

function [th,rms] = Dudkofun(fs,dF,Tmean,Fdot,theta0)
% Dudko function for puling trace
  par.nu = 1/2;
  par.model = 'DHS';
  [pd,edges] = probdens(fs,dF);
  [th,rms] = fit_Dudko_unfold(pd,edges,Tmean,Fdot,theta0,par);
end

