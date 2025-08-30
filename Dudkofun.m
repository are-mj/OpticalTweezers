
function [th,rms] = Dudkofun(fs,dF,Tmean,Fdot,theta0)
% Dudko function for puling trace
  par.nu = 1/2;
  par.model = 'DHS';
  [pd,edges] = probdens(fs,dF);
  figure;
  F = midpoints(edges);
  bar(F,pd,1);
  [th,rms] = fit_Dudko_unfold(pd,edges,Tmean,Fdot,theta0,par);
  hold on
  Fplot = linspace(edges(1),edges(end));
  DG = th(1);
  dx = th(2);
  log10k0 = th(3);
  a = par.nu*dx/DG;  
  thetacalc = [DG;a;log10k0];
  h = plot(Fplot,Dudko_unfold_probability(thetacalc,Fplot,Tmean,Fdot,par),'r');
  set(h,'linewidth',2);
end
