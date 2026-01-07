function [TRIP,TZIP] = analyse_many(files,plotting,par)
% Anslyse list of files and concatenate results tables in TRIP and TZIP
% Input:
%   files: list of file names.  
%          Either:  Full path and file name, starting with 'C:\'
%          or:      Short name, such that fullfile(datafolder,file) gives
%                   complete path and file
%   plotting:  1 - plot results
%              2 - plot temperature
%              0 or absent:  Do not plot
% SEE ALSO:: collect_tables.m

if nargin < 3
  par = params;
end
if nargin < 2
  plotting = 0;
end
TRIP = [];
TZIP = [];
for i = 1:numel(files)
  % try
    [Trip,Tzip,pull,relax,t,f,x,T] = analyse_experiment(files(i),plotting,par);
    if plotting == 2
      figure; plot(t,T)
      title(files(i));
    end
  % catch ME
  %   rethrow(ME)
  %   continue  % Skip files that give error
  % end
  TRIP =[TRIP;Trip];
  TZIP = [TZIP;Tzip];
  fprintf('Rips: %4d, Zips: %4d Filename: %s\n',height(Trip),height(Tzip),files(i));
  if plotting
    drawnow;
  end
end
  