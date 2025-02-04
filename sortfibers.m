function fiberlist = sortfibers(T)
% Fibers in optical  tweezers result tables ordered by number of rips/zips
% Input: T optical twezzers result table from collect_results/analyse_many
% Example: 
  % folders = ["20231128","20231129","20250111"];
  % files = textfile_names(folders);
  % [TP,TR] = analyse_many(files);
  % fiberlist = sortfibers([TP;TR])
  
  fiber_array = arrayfun(@(x) fiber(x),T.Filename);
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


