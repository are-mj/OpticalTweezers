function [t,f,x,shortname] = read_experiment_file(filename,detrend_x) 
% Reads a text file from molecular tweezers protein streching experient
%   Either: Pre-processed file with three columns (e.g. from CleanData app)
%   or original file from Tweezers software.
% Input: filename - including full path if not in Matlab's search path
%        detrend_x - 1: detrending of x, 0: no detrending  (default: 1)
% Output:
%    t   : Time (s) 
%    f   : Force (pN)
%    x   : Trap position (mean of two coluns A_dist_y and B_dist)
%    shortname: file name without path. Form: 'xx.txt' or '<subfolder>/xX.txt'
% Depending on the file format, time is read from Time column in file or 
%  Calculated as CycleCounts*4000;

% Author: Are Mjaavatten
% Version 2023-12-02: Switched to using readtable
% Simplified version 2024-01-20
% Version 2024-03-31: skip initial lines with negative diff(t)

  if nargin < 2
    detrend_x = 1;  % Detrending of x is default
  end
  d = dir(filename);
  if isempty(d)  % File not found in current path
    error('File not found');
  end
  filename = strrep(filename,'\','/');  % Use Unix separator
  filename_slashes = regexp(filename,'\/');  % Position of '/' in files
  fn = char(filename);  % Translate to character array
  if numel(filename_slashes) < 2
    shortname = fn;  % form: 'xx.txt' or '<subfolder>/xX.txt'
  else
    shortname = fn(filename_slashes(end-1)+1:end); % '<subfolder>/xX.txt'
  end
  shortname = string(shortname);  % Convert back to string

  warning('off','MATLAB:table:ModifiedAndSavedVarnames');
  data = readtable(filename);
  warning('on','MATLAB:table:ModifiedAndSavedVarnames');

  data = rmmissing(data);  % Remove rows with NaNs

%% Brief file format
  if width(data) == 3   
    cols = data.Properties.VariableNames;
    if any(contains(cols,'t'))&&any(contains(cols,'x')) ...
        &&any(contains(cols,'f'))
      t = data.t;
      x = data.x;
      f = data.f;     
    else
      t = data.Var1;
      x = data.Var2;
      f = data.Var3;
    end
    return
  end
%% Full file format:
  timecol = contains(data.Properties.VariableNames,'time_sec_');
  if any(timecol)
    t = data.time_sec_;
  else
    cps = 4000;  % CycleCounts per second
    countscol = contains(data.Properties.VariableNames,'CycleCount');
    if any(countscol)
      t = data.CycleCount/cps;
    else
      t = [];
      f = [];
      x = [];
      return
    end    
  end

  if numel(t) < 10
    f = [];
    x = [];
    return
  end
  start = 1;
  negdt = find(diff(t(1:10))<0);
  if ~isempty(negdt)
    start = negdt+1;  % skip any high t values at start
  end
  t = t(start:end);
  f = -data.Y_force(start:end);
  xA = data.A_dist_Y(start:end);
  xB = data.B_dist_Y(start:end);
  x = mean([xA,xB],2);
  if detrend_x
    x = detrend(x); 
  end
end
