function [t,f,x,T] = remove_time_loops(t,f,x,T)
%  Occasionally, the recorded time series contains repeated episodes. 
%  This function removes such episodes from the data

  % Eliminate records containing NaN:
  bad = isnan(x) | isnan(f) | isnan(t);
  f(bad) = []; x(bad)=[]; t(bad) = [];T(bad) = [];

  if sum(diff(t)<0)== 0  % No time repeats
    return
  end
  while sum(diff(t)<0) > 0
    loopstart = find(diff(t)<0,1);  % Last index before loop
    loopend = find(t(loopstart+1:end)>=t(loopstart),1)+loopstart;
    if isempty(loopend)  % No time > s.t(loopstart) found
      loopend = numel(t);  % Remove to end of time series
    end
    bad = loopstart:loopend;  
    t(bad) = [];
    f(bad) = [];
    x(bad) = [];
    T(bad) = [];
  end

