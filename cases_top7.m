function [ripcases,zipcases,clusters] = cases_top7(TRIP,TZIP)
% Create logical arrays for combinations of 
%   Temperature
%   Pulling speed
%   clusters
% Output:
%   ripcases : Unfolding cases
%     ripcases(m).selected is a logical column vector
%       True if TRIP.Temperature is in ripcases(m).Tclass
%       and TRIP.Pullingspeed is in ripcases(m).speedclass
%     ripcases(m).text describes tha case
%   zipcases : Refolding cases

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
  clusters = clusterdefs(Trip);
  
  m = 0;
  for i = 1:4  % Temp
    if i == 2 ||  i == 3
      speeds = 2;
    else 
      speeds = 1:3;
    end  
    for j = speeds % Speed
      m = m+1;
      ripcases(m).selected = Tclass(:,i) & speedclass(:,j);
      ripcases(m).text = strcat("Unfolding, ",Ttext(i),", ",speedtext(j));
      ripcases(m).clusters = ripcases(m).selected & clusters;
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
      else ripcases
        speeds = 1:3;
      end
      for j = speeds % Speed
          m = m+1;
          zipcases(m).selected = Tclass(:,i) & speedclass(:,j);
          zipcases(m).text = strcat("Refolding, ",Ttext(i),", ",speedtext(j));
          zipcases(m).nrips = sum(ripcases(m).selected);
      end
    end
  else
    zipcases = NaN;
  end
end