function [th,rms] = BellfunZip(fs,dF,Tmean,Fdot,theta0)
% Bell function for relaxing trace
  [pd,edges] = probdens(fs,dF);
  figure;
  F = midpoints(edges);
  bar(F,pd,1);
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