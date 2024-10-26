function [TRIP,TZIP] = analyse_many(files,plotting)
% Anslyse list of files and concatenate results tables in TRIP and TZIP
% Input:
%   files: list of file names.  
%          Either:  Full path and file name, starting with 'C:\'
%          or:      Short name, such that fullfile(datafolder,file) gives
%                   complete path and file
%   plotting:  1 - plot results
%              0 or absent:  Do not plot

if nargin < 2
  plotting = 0;
end
TRIP = [];
TZIP = [];
for i = 1:numel(files)
  if plotting
    figure(i);
    clf;
  end
  [Trip,Tzip] = analyse_experiment(files(i),parameter_struct,plotting);
  TRIP =[TRIP;Trip];
  TZIP = [TZIP;Tzip];
  fprintf('Rips: %4d, Zips: %4d Filename: %s\n',height(Trip),height(Tzip),files(i));
  if plotting
    drawnow;
  end
end
  