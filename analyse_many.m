function [TP,TR] = analyse_many(files,plotting)
% Anslyse list of files and concatenate results tables in TP and TR
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
TP = [];
TR = [];
for i = 1:numel(files)
  if plotting
    figure;
  end
  [Tp,Tr] = analyse_experiment(files(i),plotting);
  TP =[TP;Tp];
  TR = [TR;Tr];
  if plotting
    drawnow;
  end
end
  