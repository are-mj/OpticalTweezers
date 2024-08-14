function r = valid_trace_part(f,sgn)
% Remove low slope parts at ends of trace
  nf = numel(f);
  r = 1:nf;    % Full index range
  f = sgn*f;

  slope = movingslope(f,max(round(nf/50),2));
  if sgn<0
    % Test for rip delayed to relaxing trace:
    [val,loc] = max(slope);
    midpoint = floor(numel(f)/2);
    maxvalid = max(slope(midpoint:end));
    if loc < midpoint && val>maxvalid*2
      % Large peak, probably late rip.  Flatten peak, not to confuse the
      % calculation of run lengths.
      slope(slope>maxvalid)=maxvalid;
    end
  end
      
  slopemeans = movmean(slope,round(nf/10));
  lowslope1 = rle(slopemeans<max(mean(slope)/3,0));
  % Find run lengths of low and high force:
  % lowslope = rle(slope<max(slope)/4);
  lowslope = rle(slope<max(mean(slope)/3,0));

  
  % Convert brief perionds of high slope (probably due to noise) to low
  brief = round(nf/25); % 4% of total
  slopetype = lowslope(1,:);
  slopetype(~slopetype & lowslope(2,:)<brief) = 1;
  lowslope(1,:) = slopetype;
  ix = cumsum(lowslope(2,:));
  % Collapse lowslope:
  x = rle(lowslope(1,:));
  ixr = cumsum(x(2,:));
  type = lowslope(1,ixr);
  if type(1) == 1 % Starts with low force
    startrange = 1:ix(ixr(1));
  else 
    startrange = [];
  end

 
  
  % % Eliminate flat part at start:
  % % Find index for first large increase in f:
  % startsteep = find(slope>10,1);
  % if isempty(startsteep) | f(ix(1))>10
  %   startrange = [];
  % else
  %   startrange = 1:startsteep;
  % end

  % Eliminate flat parts at end:
  endrange = [];
  if lowslope(2,end)  % low slope at end
    lastlow = find(~lowslope(1,:),1,'last');
    endrange = ix(lastlow):ix(end);
  end
  
  r([startrange,endrange]) = [];
end
