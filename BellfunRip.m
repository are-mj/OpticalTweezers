function [th,rms] = BellfunRip(fs,dF,Tmean,Fdot,theta0)
% Bell function for pulling trace
  [pd,edges] = probdens(fs,dF);
  [th,rms] = fit_Bell_unfold(pd,edges,Tmean,Fdot,theta0);
  
  figure;
  F = midpoints(edges);
  bar(F,pd,1);
  hold on
  Fplot = linspace(edges(1),edges(end));
  h = plot(Fplot,Bell_unfold_probability(th,Fplot,Tmean,Fdot),'r');
  set(h,'linewidth',2);
end