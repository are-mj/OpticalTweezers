function [TRIP,TZIP] = analyse_many(files,plotting,par)
% Anslyse list of files and concatenate results tables in TRIP and TZIP
% Input:
%   files: list of file names.  
%          Either:  Full path and file name, starting with 'C:\'
%          or:      Short name, such that fullfile(datafolder,file) g
  [Trip,Tzip] = analyse_experiment(files(i),plotting,par);
  TRIP =[TRIP;Trip];
  TZIP = [TZIP;Tzip];
  fprintf('Rips: %4d, Zips: %4d Filename: %s\n',height(Trip),height(Tzip),files(i));
  if plotting
    drawnow;
  end
end
  