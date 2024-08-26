function r = valid_trace_part(f,sgn)
% Remove low slope parts at ends of trace
  nf = numel(f);
  r = 1:nf;    % Full index range
  f = sgn*f;

  noise = std(f-smoothdata(f,"movmean",round(nf/10)));
  if noise > 1  % Data too noisy.  Axxept full trace
    return
  end

  slope = movingslope(f,max(round(nf/50),2));
  % Eliminate any parts with negative slope at start:
  negslope = rle(slope<0);
  ix = cumsum(negslope(2,:));
  endnegslope = max(ix(negslope(1,:)>0 & ix<numel(f)/7));
  if isempty(endnegslope)
    endnegslope = 1;
  end  
  if sgn<0
    % Test for period of high slope (of -f)  (possible late rip):
    hislope = find(slope>(max(slope)+median(slope))/2);
    % force change during high slope:
    df = f(hislope(end))-f(hislope(1));
    % If the total force change during high slope is small, the high slope 
    % is probably a delayed rip in the relaxing trace.  Cap the slope
    % here to not confuse the search for a zip:
    if df/(max(f)-min(f)) < 0.3 
      slope(hislope) = median(slope);
    end
    % Otherwise, the high slope region is probably the real relaxation

    % [val,loc] = max(slope);
    % midpoint = floor(numel(f)/2);
    % maxvalid = max(slope(midpoint:end));
    % if loc < midpoint && val>maxvalid*2
    %   % Large peak, probably late rip.  Flatten peak, not to confuse the
    %   % calculation of run lengths.
    %   slope(slope>maxvalid)=maxvalid
    % end
  end
      
  slopemeans = movmean(slope,round(nf/10));
  % Find run lengths of low and high slope:
  % lowslope = rle(slope<max(slope)/4);
  % lowslope = rle(slope<max(mean(slope)/3,0));
  % lowslope = rle(slopemeans<max(slopemeans)/3);
  lowslope = rle(slope<max(slopemeans)/4);
  % This works poorly for very noisyd dats
  
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
    startrange = 1:max(endnegslope,ix(ixr(1)));
  else 
    startrange = 1:endnegslope;
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
