function fiberlist = sortfibers(folders)
% Fibers in folders sorted by number of rips and zips
  files = textfile_names(folders);
  [TP,TR] = analyse_many(files);
  fiber_array = arrayfun(@(x) fiber(x),[TP.Filename;TR.Filename]);
  fibers =unique(fiber_array);
  counts = zeros(length(fibers),1);
  for i = 1:length(fibers)
    counts(i)= length(find(fiber_array==fibers(i)));
  end
  T = table(fibers,counts);
  fiberlist = sortrows(T,'counts','descend');
end

function f = fiber(filename)
% Returns only the fiber part of filename
  slash = regexp(filename,'\/');
  f = extractBetween(filename,1,slash+1);
end


