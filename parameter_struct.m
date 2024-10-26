function par = parameter_struct
% parameters used by analyse_experiment.m  

% Parameters for separating out individual traces:
  par.threshold = 15; % (pN) Crossing this value defines high and low force periods
  par.lim = 20; % (pN) High force periods that stays below this value are regarded as noise

% Parameters for rip:finder:
  % Minimum and maximum fractions of trace length to be used for fitting
  % straight lines bedore and after a rip/zip
  par.minfitfraction = 0.05;  
  par.maxfitfraction = 0.15;  
  par.supportlength = 5; % number of points used for the movingslope window
  par.ripsteps = 1;      % Number of timesteps from rip start to steepest force change
  par.min_fstep = 0.4;   % Minimum reduction in force at an unfoding (pN)
  par.maxrips = 1;       % Maximum number of rips/zips accepted per trace
  par.overstretch = 55;  % Maximum rip/zip force (probably overstretch)
  par.noisefactor = [3,1]; % Skip rips/zips if fstep/noise < par.noisefactor
                           % [rip,zip]   
% Extra heating table. 
% Row 1: heater setting, from digits 2 and 3 in status column.
% Row 2: Corresponding temperature rise over mean COM file value
  % par.Tlist = [0 2 4 6 8 10 14 12 16 20 24 31; ...
  % 0 3.07 6.96 10.78 13.75 16.92 20.14 22.92 25.10 28.01 30.77 34.83];
% Ad-hoc change to handle IR laser (09=
  par.Tlist = [0 2 4 6 8 9 10 14 12 16 20 24 31; ...
  0 3.07 6.96 10.78 13.75 -100 16.92 20.14 22.92 25.10 28.01 30.77 34.83];
end

