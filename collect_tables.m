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
      % fprintf("Collects tables from the Outputfolder for %s\n",filelist(i));
      fprintf('Collects from app output. Rips: %4d, Zips: %4d Filename: %s\n',height(Trip),height(Tzip),files(i));
    else
      [Trip,Tzip] = analyse_experiment(filelist(i));
      % fprintf("Calculates tables for %s\n",filelist(i));
      fprintf('Calcuates tables.         Rips: %4d, Zips: %4d Filename: %s\n',height(Trip),height(Tzip),files(i));
    end
    if ~isempty(Trip)
      TRIP = [TRIP;Trip];
    end
    if ~isempty(Tzip)
      TZIP = [TZIP;Tzip];
    end
  end