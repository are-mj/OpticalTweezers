function r = valid_trace_part(f,sgn)
% Remove flat psrts at beginning and end of trace
  nf = numel(f);
  r = 1:nf;

  noise = std(f-smoothdata(f,"movmean",round(nf/10)));
  if noise > 1  % Data too noisy.  Accept full trace
    return
  end
  firstgood = 1;
  lastgood = nf;  
  f = sgn*f;
  n =  histcounts(f,linspace(min(f),max(f)));
  warning('off','signal:findpeaks:largeMinPeakHeight');
  [~,pks] =findpeaks([0,n,0],'minpeakheight',mean(n)*4);
  pkstart = pks(pks<10)-1; % histogram peak near start
  pkend = pks(pks>90)-1;  % histogram peak 
  if isempty(pkstart)
    pkstart = 1;
  else
    % Find any excess of data points near trace start
    % i.e. flat part of trace at start
    ixstart = find(n(pkstart(end):end)<mean(n)*2,1)+pkstart;
    if ~isempty(ixstart)
      firstgood = sum(n(1:ixstart(end)));
    end
  end
  if ~isempty(pkend)
    % Find any excess of data points near trace end
    % i.e. flat part of trace at at end    
    ixend = find(n(pkstart(1):pkend(1))<mean(n)*2,1,'last')+pkstart;
    if ~isempty(ixend)
      lastgood = nf-sum(n(ixend(1):end));
    end
  end
  warning('on','signal:findpeaks:largeMinPeakHeight');
  r = firstgood:lastgood;
end
