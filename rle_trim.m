function r = rle_trim(r,cutoff,minrun)
% Removes all run lengths with values below cutoff.
% Next, consolidataes all neighboring run lengths of identical value
% Input: r:  Output of rle
%        cutoff: Eliminate all columns where r(1,:) < cutoff
%        minrun: Disregard all runs shorter that minrun

  r = [r;cumsum(r(2,:))];     % add row with cumulative indices
  r(:,r(1,:) < cutoff) = [];  % remove entries below cutoff

  % Join short runs with the following run:
  short = find(r(2,:)<minrun);
  ncols = size(r,2);
  while ~isempty(short)
    if short(end)==ncols
      r(:,ncols)= [];
      ncols = ncols-1;
      short(end) = [];
    else
      r(:,short) = [];
    end
    short = find(r(2,:)<minrun);
  end

  % Compact r: Join all neighboring runs with same value:
  rval = rle(r(1,:));
  ixrval = [0,cumsum(rval(2,:))];
  for j = numel(ixrval):-1:2
    r(:,ixrval(j-1)+1:ixrval(j)-1) = [];
  end
end



  