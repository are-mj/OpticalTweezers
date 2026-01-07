function par = params
% parameters used by analyse_new.m  

%  SETTINGS:
% Specify which rips to look for:
  par.maxrips = 1;       % Maximum number of rips/zips accepted per trace
  par.laterips = 0;      % 1: Look for rips in the relax trace. 0: skip
  par.maxpointspertrace = 1000; % Decimate time series if the mean number of 
                         % records per trace exceeds this value.
                         % Set to 0 for no decimation

% Worm-Like-Chain parameters
  % Top7:
  par.WLC_P = 0.65;   % Persistence  length (nm)
  par.WLC_L0 = 29.28;  % Contour length (nm)  

% TUNING PARAMETERS FOR ANALYSE_EXPERIMEMT:
% Parameters for separating out individual traces:
  par.threshold = 12; % (pN) Crossing this value defines high and low force periods
  par.lim = 16; % (pN) High force periods that stay below this value are disregarded as noise

% Parameters for rip_finder:
  % Minimum and maximum fractions of trace length to be used for fitting
  % straight lines before and after a rip/zip
  par.minfitfraction = 0.05;  
  par.maxfitfraction = 0.15;  
  par.supportlength = 5; % number of points used for the movingslope window
  par.minpeak_slope = 0.003; % Minimum slope peak height in singlerip_finder
  par.ripsteps = 1;      % Number of timesteps from rip start to steepest force change
  par.min_fstep = 0.4;   % Minimum reduction in force at an unfoding (pN)
  par.overstretch = 55;  % Maximum rip/zip force (probably overstretch)
  par.noisefactor = [2,1]; % Skip rips/zips if fstep/noise < par.noisefactor
                           % [rip,zip]   
  par.maxzipfactor = 0.65; % Discard zips that occur at a higher force than
                           % par.maxzipfactor*max(force for the current
                           % trace)
% deltax limits (discard rips/zips outside limits)
  par.deltaxlimits_rips = [5,30];
  par.deltaxlimits_zips = [-30,-5];
  
% Extra heating table;  
% Row 1: heater setting, from digits 2 and 3 in status column.
% Row 2: Corresponding temperature rise over mean COM file value
% Depends on instrument name
  par.Tlist{1} = [0 2 4 6 8 10 12 14 16 20 24 31; ...
  0 3.07 6.96 10.78 13.75 16.92 20.14 22.92 25.10 28.01 30.77 34.83];
  par.Tlist{2} = [0 2 3 4 6 16;0 5.48 7.72 10.49 14.49 27.30];
  par.Tlist{3} = [0;0];
  par.Instrumentname = ["Tim's Gift 845-845 nm","SBS-tester 850-808", ...
    "BlueMini"];
  % These are the instruments in Christian's and Steve's labolratories,
  % respectively

% Ad-hoc change to handle IR laser
  % par.Tlist = [0 2 4 6 8 9 10 14 12 16 20 24 31; ...
  % 0 3.07 6.96 10.78 13.75 -100 16.92 20.14 22.92 25.10 28.01 30.77 34.83];
end

