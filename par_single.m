function par = par_single
% parameters for single rip/zip per trace:

% Parameters for separating out individual traces:
  par.min_trace_distance = 1; % Minimum distance between traces (seconds)


% Parameters for rip:finder:
  par.nus = 200;         % Sampling frequency (samples per second)

  % Minimum and maximum fractions of trace length to be used for fitting
  % straight lines bedore and after a rip/zip
  par.minfitfraction = 0.05;  
  par.maxfitfraction = 0.15;  

  par.supportlength = 5; % number of points used for the movingslope window
  par.ripsteps = 1;      % Number of timesteps from rip start to steepest force change
  par.min_fstep = 0.4;   % Minimum reduction in force at an unfoding (pN)
  par.maxrips = 1;
  par.overstretch = 55;  % Maximum rip/zip force (probably overstretch)
  par.noisefactor = [3,1]; % Skip rips/zips if fstep/noise < par.noisefactor
                           % [rip,zip]                           
end

