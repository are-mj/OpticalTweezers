function test_fitting
% Testing the fitting functions

  load Tables TRIP
  ucases = cases(TRIP);
  selection = ucases(2).clusters(:,2);  % Select data for one cluster
  force = TRIP.Force(selection);  
  dF = .333;
  [pd_obs,edges] = probdens(force,dF);
  Tmean = mean(TRIP.Temperature(selection));
  Fdot = mean(TRIP.Fdot(selection));
  theta0 = [1;-3];
  [theta,resnorm] = fit_Bell_unfold(pd_obs,edges,Tmean,Fdot,theta0);
  F = midpoints(edges);
  Fplot = linspace(edges(1),edges(end));
  pdplot = Bell_unfold_probability(theta,Fplot,Tmean,Fdot);
  figure;
  bar(F,pd_obs,1);
  hold on;
  plot(Fplot,pdplot,'r','linewidth',2);
end

function [ucases,rcases] = cases(TRIP,TZIP)
% Create logical arrays for combinations of 
%   Temperature
%   Pulling speed
%   Cluster
% Output:
%   ucases : Unfoldig cases
%     ucases(m).selected is a logical column vector
%       True if TRIP.Temperature is in ucases(m).Tclass
%       and TRIP.Pullingspeed is in ucases(m).speedclass
%     ucases(m).clusters is a n by 3 (or 2) logical array
%     ucases(m).text describes tha case
%   rcases : Refolding cases
%     Similar to ucases but lacking clusters

  if nargin < 1
    load Tables TRIP TZIP
  end
  
  % Unfold cases:
  T = TRIP.Temperature;
  Tclass = [3<T&T<=7 , 7<T&T<=14 , 14<T&T<=20,20<T&T<=30];
  Ttext = ["3-7","7-14","14-20","20-30"];
  
  speed = TRIP.Pullingspeed;
  fast = speed>250;
  normal = speed<250 & speed>50;
  slow = speed<50;
  speedclass = [slow,normal,fast];
  speedtext = ["<50","50-250",">250"];
  
  [cl1,cl2,cl3] = clusterdefinitions(TRIP);
  Clusters = [cl1,cl2,cl3]; 
  
  m = 0;
  for i = 1:4  % Temp
    if i == 2 ||  i == 3
      speeds = 2;
    else 
      speeds = 1:3;
    end  
    for j = speeds % Speed
      m = m+1;
      ucases(m).selected = Tclass(:,i) & speedclass(:,j);
      ucases(m).text = strcat("Unfolding, ",Ttext(i),", ",speedtext(j));
      ucases(m).clusters = Clusters & ucases(m).selected;
      if i==4 && j == 3
        ucases(m).clusters = [cl1,cl2|cl3] & ucases(m).selected;
      end
      ucases(m).nrips = sum(ucases(m).clusters);
    end
  end

  if nargin > 1
    % refold cases:
    T = TZIP.Temperature;
    Tclass = [3<T&T<=7 , 7<T&T<=14 , 14<T&T<=20,20<T&T<=30];
    Ttext = ["3-7","7-14","14-20","20-30"];
    
    speed = TZIP.Pullingspeed;
    fast = speed>250;
    normal = speed<250 & speed>50;
    slow = speed<50;
    speedclass = [slow,normal,fast];
    speedtext = ["<50","50-250",">250"];
    
    m = 0;
    for i = 1:4  % Temp
      if i == 2 ||  i == 3
        speeds = 2;
      else 
        speeds = 1:3;
      end
      for j = speeds % Speed
          m = m+1;
          rcases(m).selected = Tclass(:,i) & speedclass(:,j);
          rcases(m).text = strcat("Refolding, ",Ttext(i),", ",speedtext(j));
          rcases(m).nrips = sum(ucases(m).selected);
      end
    end
  else
    rcases = NaN;
  end
end

function [cl1,cl2,cl3,outliers,clustershapes] = clusterdefinitions(Trip)
% Defining clusters for Top7  
% Input: Trip: Results table for rips
% Output: Clusters: 3 by height(Trip) logical array
%         Clusters(i,j) = true if rip no. i is in cluster j
%         clustersahpes: polyshape objects for plotting the cluster extent
%         in a Force  vs Delta x scatter plot
% NOTE:  If very few rips (e.g for slow pulling speed), clusters 2 and 3
%        is often combined to 1
%        Rips that are not in any cluster are denoted outliers:
%          outliers = ~any(clusters,2)
% The clusters are defined as rectangles in a Δx - Force scatter plot
% The borders of the ractangles are adjusted so that the Bell model
% root-mean-square deviation between observation and fitted model is close
% to a minimum.

  forcegrid = [6,11,16,55];
  dxgrid = [10,20,26];
  cl1 = Trip.Deltax > dxgrid(2) & Trip.Deltax < dxgrid(3) & ...
        Trip.Force > forcegrid(3) & Trip.Force < forcegrid(4);
  cl2 = Trip.Deltax > dxgrid(1) & Trip.Deltax < dxgrid(2) & ...
        Trip.Force > forcegrid(2) & Trip.Force < forcegrid(3);
  cl3 = Trip.Deltax > dxgrid(1) & Trip.Deltax < dxgrid(2) & ...
        Trip.Force > forcegrid(1) & Trip.Force < forcegrid(2);
  outliers = ~(cl1|cl2|cl3);

  clustershapes(1) = polyshape(dxgrid([2,3,3,2,2]),forcegrid([3,3,4,4,3]));
  clustershapes(2) = polyshape(dxgrid([1,2,2,1,1]),forcegrid([2,2,3,3,2]));
  clustershapes(3) = polyshape(dxgrid([1,2,2,1,1]),forcegrid([1,1,2,2,1]));
end