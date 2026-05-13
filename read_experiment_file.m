function data = read_experiment_file(file,Tlist,detrend_x) 
% Reads a test file from Steven Smiths's minitweezer instrument
%   Can also read a file with only three columns (t,x,f)
% Input: filename - including full path if not in Matlab's search path
%        Tlist: Table specifying extra heating from the status column
%           Skip if nargout < 4
%           Tlist = []:  Read Tlist from params.m if possible
%                        Otherwise use T from COM file
%           Tlist = single number: Set temperature to this number
%           Tlist = NaN:  Report temperatures as NaN
%        detrend_x:  1: detrending of x, 0: no detrending  (default: 0)
% Output:
%    data = [t,f,x] (or [t,f,x,T] if temperature data are avaailable
%    t   : Time (s) 
%    f   : Force (pN)
%    x   : Trap position (mean of two coluns A_dist_y and B_dist)
%    T   : Temperature (°C)
%
% Depending on the file format, time is read from Time column in file or 
%  Calculated as CycleCounts*4000;

% Author: Are Mjaavatten

  data = [];  % make sure the output is defined

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
  end
  if ~isfile(filename)
    error("File %s is not found",filename);
  end
  filename = strrep(filename,'\','/');  % Use Unix separator
  warning('off','MATLAB:table:ModifiedAndSavedVarnames');
  indata = readtable(filename);
  warning('on','MATLAB:table:ModifiedAndSavedVarnames');

  indata = rmmissing(indata);  % Remove rows with NaNs or missing values

%% Brief file format
  if width(indata) == 3   
    cols = indata.Properties.VariableNames;
    if any(contains(cols,'t'))&&any(contains(cols,'x')) ...
        &&any(contains(cols,'f'))
      t = indata.t;
      x = indata.x;
      f = indata.f;     
    else
      t = indata.Var1;
      x = indata.Var2;
      f = indata.Var3;
    end
    data = [t,f,x];
    data(any(isnan(data),2),:) = [];
    return
  end
%% Full file format (Mini Tweeezers from Steven B. Smith)
  timecol = contains(indata.Properties.VariableNames,'time_sec_');
  if any(timecol)
    t = indata.time_sec_;
  else
    cps = 4000;  % CycleCounts per second
    countscol = find(contains(indata.Properties.VariableNames,'CycleCount'));
    if any(countscol)
      % t = indata.CycleCount/cps;
      t = table2array(indata(:,countscol))/cps;
    end    
  end

  startrows = 1:min(10,length(t));
  negdt = find(diff(t(startrows))<0);
  if ~isempty(negdt)
    start = negdt+1;  % skip any high t values at start
  else
    start = 2;  % Skip first record, which seldom makes sense
  end
  t = t(start:end);
  f = -indata.Y_force(start:end);
  xA = indata.A_dist_Y(start:end);
  xB = indata.B_dist_Y(start:end);
  status = indata.Status(start:end);  
  x = mean([xA,xB],2);
  if numel(t) < 10
    T = NaN(size(t));
    data = [t,x,f,T];
    % warning('Too few rows in %s\n',filename)
    return
  end
  if detrend_x
    x = detrend(x); 
  end
  data = [t,f,x];
  data(any(isnan(data),2),:) = [];

  % *** Temperature  ***  
  T_OK = false;
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
        instrumentno = find(strcmp(instrument,par.Instrumentname));
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
  data = [data,T];
end
