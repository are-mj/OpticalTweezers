function files = textfile_names(folders)
% Finds all *.txt files in the string array of folder names
% skips file names containing 'COM'
  if ischar(folders)
    folders = string(folders);
  end
  nfolders = numel(folders);
  files = [];
  for i = 1:nfolders
    d = dir(fullfile(datafolder,folders(i),'*.txt'));
    nfiles = numel(d);
    newfiles = [];
    for j = 1:nfiles
      name = string(d(j).name);
      if ~contains(name,"COM")
        newfiles =[newfiles;strcat(folders(i),"/",name)];
      end
    end
    files = [files;newfiles];
  end
end