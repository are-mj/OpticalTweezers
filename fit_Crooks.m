function [DG,DGci,Wsu_kcal,Wu_kcal,Wsr_kcal,Wr_kcal] = fit_Crooks(TRIP,TZIP,par,plotting)
% Calculate DG from unfolding and refolding parameters dx and log10k0 
% DG equals the work at the point where unfolding and refolding work are
% equal, according to the Crooks fluctuation theorem
% Input:  
%   TRIP, TZIP: Output tables from analyse_experiment or analyse_many
%   plotting: show plots. Default: no plots
% Output:
%   DG:    Crooks DeltaG
%   DGci:  Confidence interval calculated by bootstrapping

% 20240209: Specified analytical search range in match.
% 20250529: Changed from standard deviation to 95% confidence interval
% 20250724: Modified to handle a single data set only

  conversion = 0.1439326;  % Energy units kcal/kmol  

  if nargin < 4
    plotting = 0;
  end

  L0 = par.WLC_L0;
  P = par.WLC_P;
  Tmean = mean(TRIP.Temperature,"omitmissing");
  T = Tmean + 273.15;  % °C to K

  % Refold:
  f_refold = TZIP.Force;
  bad = f_refold < 2;
  f_refold(bad) = [];  % Skip unrealistic values
  deltax_refold = -TZIP.Deltax;
  deltax_refold(bad) = []; 
  Wsr = stretchwork(f_refold,deltax_refold,P,T,L0); % Should we subtract Wsr here?
  Wr = f_refold.*deltax_refold - Wsr;
  % Convert to kcal/mol to get the correct pdr plot
  Wr_kcal = Wr*conversion;  % Convert energy units
  Wsr_kcal = Wsr*conversion;
  pdr = fitdist(Wr_kcal,'normal');   % Normal distribution

  % Unfold:
  f_unfold = TRIP.Force;
  deltax_unfold = TRIP.Deltax;

  if plotting
    h = plot(pdr);
    h(1).Color = 'k';
    hold on;
    % colors = [0 0 1;0.5 0.2 0.5;0.3 0.6 0.6;.5*[1 1 1]];
    % h(1).Color = colors(1,:);
  end

  fu = f_unfold;
  dxu = deltax_unfold;
  Wsu = stretchwork(fu,dxu,P,T,L0);
  Wu = fu.*dxu-Wsu;        % Net work
  Wu_kcal = Wu*conversion;
  Wsu_kcal = Wsu*conversion;
  pdu = fitdist(Wu_kcal,'normal');
  DG = match(pdr,pdu);  % kcal/mol
  if nargout > 1
    % Perform this time-consuming calculation only if needed:
    DGci = Crooks_ci(pdu,pdr);
  end
  if plotting
    h = plot(pdu);
    h(1).Color = 'k';
    % h(1).Color = colors(1,:);
    xlabel('Work (kcal/mol)');ylabel('Probalility density'); 
    legend('Refold','','Unfold','Fitted normal distributions')
  end
end

function pd = pdfun(pdobj,x)
  mu = pdobj.mu;
  sigma = pdobj.sigma;
  pd = exp(-0.5*((x-mu)/sigma).^2)/sigma/sqrt(2*pi);;
end

function DG = match(pd1,pd2)
% Find the point where two normal distributions are equal
  fun = @(x) pdfun(pd1,x)-pdfun(pd2,x);
  try
    [DG,~,exitflag] = fzero(fun,[pd1.mu,pd2.mu]);
  catch  % No crossong point between distribution tops
    [DG,~,exitflag] = fzero(fun,[pd1.mu-pd1.sigma,pd2.mu]);
  end
  if exitflag < 0
    warning('fzero problem')
  end
end

function ci = Crooks_ci(pdu,pdr)
% Monte Carlo calculation of the confidence interval of the intersection of pdr and pdu
  
  pdu_ci = pdu.paramci;
  pdu_mu_std = pdu.mu - pdu_ci(1,1);
  pdu_sigma_std = pdu.sigma - pdu_ci(1,2);

  pdr_ci = pdr.paramci;
  pdr_mu_std = pdr.mu - pdr_ci(1,1);
  pdr_sigma_std = pdr.sigma - pdr_ci(1,2);

  nsamples = 1000; 

  DG = zeros(nsamples,1);
  for i = 1:nsamples
    try  % normrnd may occasionally return negative values
         % This will crash makedist
      umu = normrnd(pdu.mu,pdu_mu_std);
      usigma = normrnd(pdu.sigma,pdu_sigma_std);
      rand_pdu = makedist('normal','mu',umu,'sigma',usigma);
      % plot(rand_pdu);
  
      rmu = normrnd(pdr.mu,pdr_mu_std);
      rsigma = normrnd(pdr.sigma,pdr_sigma_std);
      rand_pdr = makedist('normal','mu',rmu,'sigma',rsigma); 
      % plot(rand_pdr);
      DG(i) = match(rand_pdr,rand_pdu);
    catch
      DG(i) = NaN;
    end
  end
  DG(isnan(DG)) = [];
  pd = fitdist(DG,'normal');
  ci = norminv([0.025 0.975],pd.mu,pd.sigma);
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