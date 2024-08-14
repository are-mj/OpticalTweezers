function [splits,groups] = frequency_split(f)
% Finds indices where the oscillation freuency of noisy signal f changes 
% by an order of magnitute
  threshold = mean(f);
  % Make sure no f value == threshold:
  exact_hit = f==threshold;
  f(exact_hit) = f(exact_hit)*(1+5*eps); 

  % Find all points where f crosses the thresold value:
  up = find(f(2:end)-threshold>0 & f(1:end-1)-threshold<0);
  down = find(f(2:end)-threshold<0 & f(1:end-1)-threshold>0);
  
  nup = numel(up);
  ndown = numel(down);
  upmax = min(nup,ndown);
  downmax = min(ndown,nup);

  % Attempting to find run lengths of widely differing pulling frequency
  % Distance between crossings:
  g = up(1:upmax)-down(1:downmax);
  lg = log10(g);
  window = 100;
  h = movmax(lg,window);  % Moving maximum of up
  % Assign max up - down distances to bins around 100, 1000 and 10000
  % recordings:
  groups = [100,1000,10000];  % up-down distances
  lggroups = log10(groups);
  edges =[lggroups-0.5,lggroups(end)+0.5];
  [~,~,bins] = histcounts(h,edges); 
  rlebins = rle(bins);
  % Group very small runs with next run:
  small = find(rlebins(2,:) < numel(h)/100);   % negligible run length
  if ~isempty(small) & small(end)<size(rlebins,2) 
    rlebins(2,small+1) = rlebins(2,small+1) + rlebins(2,small);
    rlebins(:,small) = [];
  end
  levels = groups(rlebins(1,:));
  ix = [1,cumsum(rlebins(2,:))];
  nsplits = numel(levels)-1;
  splits = zeros(1,nsplits+2);
  splits(1)=1;
  splits(end) =numel(f);
  for i = 1:nsplits
    if levels(i) < levels(i+1) % Find first up value > sqrt(10)*levels(i)
      % first significant up
      up1 = up(ix(i+1)+find(g(ix(i+1):ix(i+2))>levels(i)*sqrt(10),1)-1);
      up2 = up(find(up>up1+levels(i+1)/2,1));
    else
      % last significant up before frequency change
      up1 = up(ix(i)+find(g(ix(i):ix(i+1))>levels(i)/sqrt(10),1,'last')-1);
      % next significant up:
      up2 = up(find(up>up1+levels(i)/2,1));
    end
    splits(i+1) = round(mean([up1,up2]));
  end
end