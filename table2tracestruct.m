function s = table2tracestruct(app,T,row,sgn)
% Create a trace struct for the current file from a results table row
  s.file = fullfile(app.Datafolder,app.Datafile);
  shortname = shorten_filename(app.file);
  if ~isequal(T.Filename(row),shortname)
    warning('Results file does not match the current file')
  end

  ix = find(app.filedata.t == T.time);
  if sgn > 0
    trace_start = app.filedata.valleys(find(app.filedata.valleys<ix,1,'last'));
    trace_stop = app.filedata.peaks(find(app.filedata.peaks>ix,1));
  else
    trace_start = app.filedata.peaks(find(app.filedata.peaks<ix,1,'last'));
    trace_stop = app.filedata.valleys(find(app.filedata.valleys>ix,1));
  end
  range = valid_trace_part(f(trace_start:trace_stop,sgn));
  s.t = app.filedata.t(range);
  s.f = app.filedata.f(range);
  s.x = app.filedata.x(range);
  s.T = app.filedata.T(range);

  s.ripx = T.Trapx(row);
  s.force = T.Force(row);
  s.fstep = T.Forceshift(row);
  s.deltax = T.Deltax(row);
  s.time = T.Time(row);
  s.fdot = T.Fdot(row);
  s.slope = T.slope(row);  
  s.rip_index = ix - range(1) +1;
  s.dt = T.Timestep(row);
  s.temperature = T.Temperature(row);
  s.noise = T.Noise(row); 
  s.pullingspeed = T.Pullingspeed(row);

  % fitting range before rip/zip
  nf = numel(range);
  fitb_start = max(ix - round(app.par.maxfitfraction*nf),range(1));
  fitrange_b = fitb_start:ix;
  s.pfx_b = polyfit(x(fitrange_b),f(fitrange_b),1);  
  % fitting range after rip/zip
  fita_stop = min(ix+2*app.par.ripsteps+round(app.par.maxfitfraction*nf),range(end));
  fitrange_a = (ix + +2*app.par.ripsteps):fita_stop;
  s.pfx_a = polyfit(x(fitrange_a),f(fitrange_a),1);
end
  

  