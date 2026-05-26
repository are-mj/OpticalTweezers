function r = valid_trace_part(f,sgn,par)
% Version 20260316.  Better at handling rips near peak x 
% Remove flat parts at beginning of pull trace and end of relax trace
  
  nf = numel(f);
  r = (1:nf)';

  % First remove any points with f > par.overstretch
  if sgn > 0
    highf = find(f>par.overstretch,1);
    if ~isempty(highf)
      f(highf:end) = [];
      r = 1:length(f);
    end
  else
    highf = find(f>par.overstretch,1,'last');
    if ~isempty(highf)
      f(1:highf) = [];
      r = 1:length(f);
    end
  end
  nf = length(f);
  
  f = sgn*f;
  [n,~,bin] =  histcounts(f,100);  % 100 bins
  nn = [0 n 0];

  [pks,loc,w] =findpeaks([0,n,0],'MinPeakProminence',100);
  w = ceil(w);

  if isempty(pks)
    return
  end
  if length(pks) > 2  % Too many flat parts. Discard trace
    r = [];
    return
  end
  % if any(loc >= 8 & loc < 95)
  %   r = [];  % Flat part not near start or end
  %   return
  % end
  remove = [];
  for pkno = 1:length(loc)
    if loc(pkno)>95
      flatpoints = find(bin>=loc(pkno)-w(pkno));
    else
      flatpoints = find(bin<=loc(pkno)+w(pkno));
    end
    % make sure flatpoints includes end points if relevant:
    if nf-flatpoints(end) < nf/20
      flatpoints(end)  = nf;
    end
    if flatpoints(1) < nf/20
      flatpoints(1) = 1;
    end
    bad = min(flatpoints):max(flatpoints);
    remove = union(remove,bad);
    if loc >= 95
      break   % merge all locs > 95
    end
  end
  r(remove) = [];
end

