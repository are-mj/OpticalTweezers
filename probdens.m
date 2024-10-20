function [pd,F,n] = probdens(force,dF)
% Histogram of probability densities from experimennts
% Similar to probability_density, but:
%   - uses force vector as input rather that Matlab Table Tun or Tre
%   - does not return Tmean or Fplot
% Input:
%  force:  Array of recorded transition forces
%  dF:   Force value spacing
% Output: 
%  pd: Column vector of probality desities.  The probability of unfolding 
%      (or refolding) in the interval <F(i)-dF/2,F(i+dF/2)> is pd(i)*dF
%  n: nUmber of rows in force
%  F: Force vector (useful if F is not specified in the call)

  if isempty(force)  % Empty cluster
    pd = NaN;
    F = NaN;
    n = 0;
    return
  end
  edges = floor(min(force)-0.5*dF)+0.5*dF:dF:ceil(max(force) + 0.5*dF)-0.5*dF;
  F = (edges(1:end-1)+edges(2:end))'/2;
  Values = histcounts(force,edges);
  n = sum(Values);
  pd = Values'/n/dF;  % Probability density
end
