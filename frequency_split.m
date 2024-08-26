function [splits,groups] = frequency_split(t,f)
% Finds indices whwere the expeiment time series should be split into 
% sections of the same oscillation frequenciy (within order of magnitude) 
% Also splits on order of magnintude changes in sampling frequency 

  % First split on sampling frequency:
  dt = diff(t);
  dtclass = round(log10(dt));
  rle_dt = rle(dtclass);
  cutoff = -2;  % disregard sampling intervals shorter than 3 ms. 
  % Handle cases with extremely short sampling intervals
  while sum((dtclass>cutoff)/numel(dt))<0.1
    cutoff = cutoff -1;
  end
  rle_dt_trimmed = rle_trim(rle_dt,cutoff,5);
  tsplits = rle_dt_trimmed(3,:)+1;

  % Next split on oscillation frequency
  threshold = mean(f);
  % Make sure no f value == threshold:
  exact_hit = f==threshold;
  f(exact_hit) = f(exact_hit)*(1+5*eps); 

  % Find all points where f crosses the thresold value:
  up = find(f(2:end)-threshold>0 & f(1:end-1)-threshold<0);
  down = find(f(2:end)-threshold<0 & f(1:end-1)-threshold>0);

  upfirst = up(1)<down(1);
  
  nup = numel(up);
  ndown = numel(down);
  upmax = min(nup,ndown);
  downmax = min(ndown,nup);

  % Attempting to find run lengths of widely differing pulling frequency
  % Distance between crossings:
  if upfirst
    g = down(1:downmax)-up(1:downmax);
  else
    g = up(1:upmax)-down(1:upmax);
  end
  lg = log10(g);
  window = 100;
  h = movmax(lg,window);  % Moving maximum of up
  % Assign max up - down distances to bins around 100, 1000 and 10000
  % recordings:
  % groups = [10,100,1000,10000];  % up-down distances
  % lggroups = log10(groups);
  edges = 1:4;   % i.e. [10,100,1000,10000]
  lggroups = (edges(2:end)+edges(1:end-1))/2; 
  groups = 10.^lggroups;
  [~,~,bins] = histcounts(h,edges); 
  rlebins = rle(bins);
  % Group very small runs with next run:
  small = find(rlebins(2,:) < numel(h)/100);   % negligible run length
  if ~isempty(small) & small(end)<size(rlebins,2) 
    rlebins(2,small+1) = rlebins(2,small+1) + rlebins(2,small);
    rlebins(:,small) = [];
  end
  levels = groups(rlebins(1,rlebins(1,:)>0));
  ix = [1,cumsum(rlebins(2,:))];
  nsplits = numel(levels)-1;
  splits = zeros(1,nsplits+2);
  splits(1)=1;
  splits(end) =numel(f);
  for i = 1:nsplits
    if levels(i) < levels(i+1) % Find first up value > sqrt(10)*levels(i)
      % find first significant up after step up in h
      up1 = up(ix(i+1)+find(g(ix(i+1):ix(i+2))>levels(i)*sqrt(10),1)-1);
      % up2 = up(find(up>up1+levels(i+1)/2,1));
    else
      % last significant up before frequency change
      up1 = up(ix(i)+find(g(ix(i):ix(i+1))>levels(i)/sqrt(10),1,'last')-1);
      % next significant up:
      % up2 = up(find(up>up1+levels(i)/2,1));
    end
    % splits(i+1) = round(mean([up1,up2]));
    splits(i+1) = up1;
  end
  splits = union(splits,tsplits);
end