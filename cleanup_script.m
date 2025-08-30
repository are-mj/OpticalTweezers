for i = 1:length(files)
  matfile = fullfile(datafolder,strrep(files(i),".txt",".mat"));
  if exist(matfile,'file')
    [~,~,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_experiment(files(i));
    load(matfile,'Trip','Tzip','limits');
    if ~isempty(Trip)
      Trip = Trip(Trip.Time > t(limits(1)) & ...
          Trip.Time < t(limits(2)),:);
    end
    if ~isempty(Tzip)
      Tzip = Tzip(Tzip.Time > t(limits(1)) & ...
        Tzip.Time < t(limits(2)),:);
    end
    save(matfile,'Trip','Tzip','limits');
  end
end
  %   listOfVariables = who('-file',matfile);
  %   if ~any(ismember(listOfVariables,'limits'))
  %     disp(i)
  %   end
  % end
% end

%   cleanup_mat(fullfile(datafolder,matfile));
% end
% function cleanup_mat(matfile)
%   if exist(matfile,"file")
%     load(matfile,'Trip','Tzip');
%     if any("Fitrange" == string(Trip.Properties.VariableNames))
%       return
%     end
%     ripnos = height(Trip);
%     if ripnos > 0
%       Trip.Filename = repmat(shorten_filename(Trip.Filename(1)),[ripnos,1]);
%       newripcols = zeros(ripnos,2);
%       Trip = addvars(Trip,newripcols,'After',"Noise",'NewVariableNames',"Fitrange");
%     end
%     zipnos = height(Tzip);
%     if zipnos > 0
%       Tzip.Filename = repmat(shorten_filename(Tzip.Filename(1)),[zipnos,1]);
%       newzipcols = zeros(height(Tzip));
%       Tzip = addvars(Tzip,newzipcols,'After',"Noise",'NewVariableNames',"Fitrange");    
%     end
%     save(matfile,"Trip","Tzip");
%   end
% end