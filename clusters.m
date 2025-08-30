function [Cluster1,Cluster2,outliers,par] = clusters(TP,par,plotting)
% Separate data points in a Deltax,Force plot into clusters and outliers
% Cluster1: logical array, true for Cluster1 rows in TP
% Cluster2: logical array, true for Cluster2 rows in TP 
% outliers: true for points outside clusters
% Input:
%   TP  Matlab table from analyse_eperiment or analyse_many
%   par: Parameters for function dbscan (Statistics and machine 
%        learning toolbox) plus scaling factor for Deltax values
%        See dbscan for details.

  if nargin < 2 || isempty(par)
    par.epsilon = 5.5;  % Neighbour distance
    par.minpts = 100;   % Minimum number of points in cluster
    par.scaling = 3;
  end


  X = [TP.Deltax*par.scaling,TP.Force];  % Scaled
  labels = dbscan(X,par.epsilon,par.minpts);
  Cluster1 = labels == 1;
  Cluster2 = labels == 2;
  outliers = labels < 1;
  if nargin>2
    if plotting
      % figure;
      plot(TP.Deltax(Cluster1),TP.Force(Cluster1),'.b', ...
        TP.Deltax(Cluster2),TP.Force(Cluster2),'.r');
      hold on;
      plot(TP.Deltax(outliers),TP.Force(outliers),'.','color',0.65*[1 1 1])
      hold off
    end
  end
end