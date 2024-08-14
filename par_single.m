function par = par_single
% parameters for single rip/zip per trace:

% Parameters for separating out individual traces:
  par.min_trace_distance = 1; % Minimum distance between traces (seconds)


% Parameters for rip:finder:
  par.maxrips = 0.1;

  par.nus = 200;         % Sampling frequency (samples per second)
  par.smoothfraction = 0.1; % smoothing window as fraction of trace length

  % Minimum and maximum fractions of trace length to be used for fitting
  % straight lines bedore and after a rip/zip
  par.minfitfraction = 0.07;  
  par.maxfitfraction = 0.1;  

  par.supportlength = 5; % number of points used for the movingslope window
  par.ripsteps = 1;      % Number of timesteps from rip start to steepest force change
  par.min_fstep = 0.4;   % Minimum reduction in force at an unfoding (pN)
  par.maxrips = 1;
  par.linespan = 50;     % Maximum point span for plotted slope lines
  par.overstretch = 55;  % Maximum rip/zip force (probably overstretch)
  par.cycleforce = 10;   % Level of smoothed force to identify cycles (pN)
  par.noisefactor = 1.5; % Skip rips/zips if fstep/noise < par.noisefac

  % Parameters for relaxing trace:
  par.relax.stdwindow = 300; 
  % par.relax.min_peak = 0.02;
  % par.relax.supportlength = 15;  
end

