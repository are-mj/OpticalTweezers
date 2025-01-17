function [t,f,x,T] = read_experiment_file(file,Tlist,detrend_x) 
% Reads a test file from Steven Smiths's minitweezer instrument
%   Can also read a file with only three columns (t,x,f)
% Input: filename - including full path if not in Matlab's search path
%        Tlist: Table specifying extra heating from the status column
%           Tlist = []:  Read Tlist from params.m if possible
%                        Otherwise use T from COM file
%           Tlist = single number: Set temperature to this number
%           Tlist = NaN:  Report temperatures as NaN
%        detrend_x:  1: detrending of x, 0: no detrending  (default: 0)
% Output:
%    t   : Time (s) 
%    f   : Force (pN)
%    x   : Trap position (mean of two coluns A_dist_y and B_dist)
%    T   : Temperature (°C)
%
% Depending on the file format, time is read from Time column in file or 
%  Calculated as CycleCounts*4000;

% Author: Are Mjaavatten

  % Make sure all outputs are defined:
  t = [];
  f = [];
  x = [];
  T = [];

  if nargin < 3
    % detrend_x = 1;  % Detrending of x is default
    detrend_x = 0;  % No detrending of x is default
  end
  if nargin < 2
    Tlist = [];
  end
  % Allow file name containing full path
  if isfile(file)
    filename = file;
  else
    filename = fullfile(datafolder,file);
    if ~isfile(filename)
      error("File %p is not found",filename);
    end
  end  
  filename = strrep(filename,'\','/');  % Use Unix separator
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
    end    
  end

  if numel(t) < 10
    return
  end
  start = 2;  % Skip first record, which seldom makes sense
  % start = 1;
  negdt = find(diff(t(1:10))<0);
  if ~isempty(negdt)
    start = negdt+1;  % skip any high t values at start
  end
  t = t(start:end);
  f = -data.Y_force(start:end);
  xA = data.A_dist_Y(start:end);
  xB = data.B_dist_Y(start:end);
  status = data.Status(start:end);  
  x = mean([xA,xB],2);
  if detrend_x
    x = detrend(x); 
  end

  % *** Temperature  ***
  try
    % Read bath temperature from COM file.
    [Tbath,instrument] = T_from_COM(filename); % Temperature outside cell
    T = ones(size(t))*Tbath;
  catch
    return
  end
  if isempty(Tlist)  % Try reading from params.m
    if exist("params.m","file")
      par = params;
      if isfield(par,'Tlist') && isfield(par,"Instrumentname")
        instrumentno = find(contains(instrument,par.Instrumentname));
        if isempty(instrumentno)
          error('Unknown instrument: %s. Cannot determine temperature',instrument);
        else
          Tlist = par.Tlist{instrumentno};
        end
      end
    else
      error('Parameter function params.m not found')
    end
  end

  if isscalar(Tlist)  % This option also handles Tlist == NaN						   
    T = Tlist*ones(size(t)); % Fixed T specified
    return
  end
 
  if ~isnan(T)
    % Read heater setting from digits 2 and 3 in the status column 
    heater_setting = floor(rem(status,1000)/10);  % Number from digits 2 and 3
    heater_setting(status<1000) = NaN; 
    for ii = 1:size(Tlist,2)
      ix = heater_setting==Tlist(1,ii);
      T(ix) = T(ix) + Tlist(2,ii);				   
    end
  end
end
