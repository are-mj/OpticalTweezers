function [trace1,trace2] = adjust_x(trace1,trace2)
% shift the x values in structs for neighbouring traces trace1 and trace2
% so that they have identical end points. trace2 must be the trace 
% immediately following trace1.  trace1 may be a pulling or a relaxing trace

  x1 = trace1.x(1);
  x2 = trace1.x(end);
  z1 = trace2.x(1);
  z2 = trace2.x(end);
  sgn = sign(trace1.f(end)-trace1.f(1));
  % Normalise so that both traces span the same range of x values
  if sgn > 0  % Pulling followed by relaxing
    xspan = (abs(x1-x2)+abs(z1-z2))/2;  
    shift1 = -x1;
    scale1 = xspan/abs(x2-x1);  
    trace1 = shiftandscale(trace1,shift1,scale1);
    shift2 = -z2;
    scale2 = xspan/abs(z1-z2);
    trace2 = shiftandscale(trace2,shift2,scale2);
  else
    xspan = (abs(x1-x2)+abs(z1-z2))/2;  
    scale1 = xspan/abs(x1-x2);
    shift1 = -x2;  
    trace1 = shiftandscale(trace1,shift1,scale1);
    shift2 = -z1;
    scale2 = xspan/abs(z1-z2);
    trace2 = shiftandscale(trace2,shift2,scale2);    
  end
end

function trace = shiftandscale(trace,shift,scale)
  trace.x = (trace.x+shift)*scale;
  if isfield(trace,'ripx')
    trace.ripx = (trace.ripx+shift)*scale;
    trace.pfx_b = shiftpoly(trace.pfx_b,shift,scale);
    trace.pfx_a = shiftpoly(trace.pfx_a,shift,scale);
  end
end

function ps = shiftpoly(p,shift,scale)
  % Shifts and scale linear polynomial p
  ps = [p(1)/scale,p(2)-p(1)*shift];
end