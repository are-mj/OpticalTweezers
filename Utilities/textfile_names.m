function files = textfile_names(folders,type)
% Finds all *.txt files in the string array of folder names
% Inout: folders: string arrsy of folder names in datafolder
%        type:    "all": all *.txt files
%                 "nocom": Skip COM files (default)

  if nargin < 2
    type = "nocom";
  end
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
      if type == "nocom" && contains(name,"COM")
        continue
      else
        newfiles =[newfiles;strcat(folders(i),"/",name)];
      end
    end   
    files = [files;newfiles];
  end
end