function [pd_obs,edges,n_obs] = probdens(force,dF)
% Histogram of probability densities from experimennts
% Input:
%  force:  Array of recorded rip or zip forces
%  dF:     Force value spacing
%  n_obs   Number of rips or zips
% Output: 
%  pd_obs: Column vector of probality densities. 
%  edges:  Force bin edges. pd_obs(i) is the probability density for unfolding 
%          (or refolding) in the interval from edges(i) to edges(i+1)
%          
%  n_obs:  Number of observed rips or zips

  if isempty(force)  % Empty cluster
    pd_obs = NaN;
    edges = NaN;
    n_obs = 0;
    return
  end
  edges = floor(min(force/dF))*dF:dF:ceil(max(force/dF))*dF;
  Values = histcounts(force,edges);
  n_obs = sum(Values);
  pd_obs = Values'/n_obs/dF;  % Probability density
end
