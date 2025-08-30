function [TRIP,TZIP] = collect_tables(filelist,Outputfolder)
% Creates total rip and zip results tables for the files in filelist
% Uses tables saved from the RipAnalysis app if they exist
% Otherwise calcuates tables using function analyse_experiment
% Input: 
%   filelist: Strnng array of file names
%   Outputfolder:  Folder used by the RipAnalysis app. 
%      Default: Uses current Outputfolder from the RipAnalysis app     
%

  if nargin < 2
    load("RipAnalysis_settings.mat","appsettings")
    Outputfolder = appsettings.Outputfolder;
  end
  TRIP = [];
  TZIP = [];
  for i = 1:numel(filelist)
    matfile = fullfile(Outputfolder,strrep(filelist(i),".txt",".mat"));
    if exist(matfile,"file")
      load(matfile,"Trip","Tzip");
      fprintf('File: %25s, Rips: %4d, Zips: %4d, Modified in app\n', filelist(i),height(Trip),height(Tzip));
    else
      [Trip,Tzip] = analyse_experiment(filelist(i));
      fprintf('File: %25s, Rips: %4d, Zips: %4d, Unmodified\n', filelist(i),height(Trip),height(Tzip));
    end
    if ~isempty(Trip)
      try
        TRIP = [TRIP;Trip];
      catch  % Modified table not compatible
        [Trip,Tzip] = analyse_experiment(filelist(i));
        fprintf('File: %25s, Rips: %4d, Zips: %4d, Modified table incompatible\n', filelist(i),height(Trip),height(Tzip));        
        TRIP = [TRIP;Trip];
      end
    end
    if ~isempty(Tzip)
      try
        TZIP = [TZIP;Tzip];
      catch
        % skip
      end
    end
  end