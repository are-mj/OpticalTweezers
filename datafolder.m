function folder = datafolder
% Returns the full path folder where the optical tweezer experiment 
% measurements are stored. 
% Example (windows)
%  folder = "C:\users\are\data\tweezers"

% This works for me on both my PCs that use diferent user names
  home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
  localfolder = '\OneDrive\Data\Chile';
  folder = [home,localfolder];
  folder = "C:\Users\are\Berkeley";
% folder ='C:\Users\mjaav_hi6z2cj\OneDrive\Data\Berkeley\B101_A7C_N59C\2010-1-25_B101\2010-1-25_B101 raw data';'
% folder = "C:\Users\mjaav_hi6z2cj\OneDrive\Documents\Tmp";
% folder = "C:\Users\are\OneDrive\Documents\Tmp";
end
