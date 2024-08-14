function folder = datafolder
% Returns the folder where the optical tweezer experiment measurements
% are stored.  
%
% Example: To plot data from file kA.txt in folder 02042022:
%  plotfile(fullfile(datafolder,'02042022','kA.txt'));
% or: 
%  plotfile(fullfile(datafolder,'02042022/kA.txt'));
  home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
  localfolder = '\OneDrive\Data\Chile';
  % localfolder = '\Dropbox\Projects\Chile\Karina_New\Clean_files';
  folder = [home,localfolder];
end
