function folder = datafolder
% Returns the full path folder where the optical tweezer experiment 
% measurements are stored. 
% Example (windows)
%  folder = "C:\users\are\data\tweezers"

% This work for me on both my PCs that use diferent user names
  home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
  localfolder = '\OneDrive\Data\Chile';
  folder = [home,localfolder];
end
