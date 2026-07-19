function st = trim_trace(data,st,par)
% Remove structs with unlikely variable values from trace struct 
% Current version checks only that deltax is within limits
  for i = length(st):-1:1
    if isempty(st(i).force)
      return
    end
    sgn = sign(st(i).fdot);
    f = data(st(i).fitrange(1):st(i).fitrange(2),2);
    bad = false;
    if sgn> 0
      deltaxlim = par.deltaxlimits_rips;
    else
      if st(i).deltax > 0  % late rip
        deltaxlim = par.deltaxlimits_rips;
      else
        deltaxlim = abs(flip(par.deltaxlimits_zips));
      end
    end
    bad = bad | abs(st(i).deltax) < deltaxlim(1) | abs(st(i).deltax) > deltaxlim(2);
    if sum(bad)> 0
      st(i) = [];
    end
  end
  if isempty(st)
    st = empty_trace;
  end
end