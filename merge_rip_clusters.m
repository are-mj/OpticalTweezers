function rip_index = merge_rip_clusters(rip_index,quality,epsilon)
% Merge rips within a cluster with minimum internal rip_index distances epsilon.
% Output: rip indices of highest-quality rip in each cluster

% Thanks to chatGPT

  x = rip_index(:);
  [xs, order] = sort(x);
  groupID_sorted = zeros(size(xs));
  gid = 1;
  groupID_sorted(1) = gid;
  
  for i = 2:length(xs)
      if xs(i) - xs(i-1) <= epsilon
          groupID_sorted(i) = gid;
      else
          gid = gid + 1;
          groupID_sorted(i) = gid;
      end
  end
  
  groupID = zeros(size(x));
  groupID(order) = groupID_sorted;

  bestPos = accumarray(groupID, (1:length(groupID))', [], ...
    @(idx) idx( quality(idx) == max(quality(idx)) ));
  rip_index = rip_index(bestPos);
end