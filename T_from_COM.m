function T = T_from_COM(file)
% Read temperature from the COM file corresponding to the experiment file
%   Backup solution if the temperature is not coded in the staus variable
%   2023-10-16: Modified the construction of COM file name.  Now works for
%   abA.txt (but no longer for aAB.txt)
%   2023-12-07: Read only TemperatureB (as recommended by Steve B. Smith)
  file = char(file);
  slashes = regexp(file,'[\/\\]');
  if numel(slashes)<2  % Short version of filename
    file_full = char(fullfile(datafolder,file));  % datafolder.m must return the folder
    file_full = strrep(file_full,'\','/');
    slashes = regexp(file_full,'[\/\\]'); 
  else
    file_full = file;
  end
  strrep(file_full,'\','/');
  
  fiber = file_full(slashes(end)+1:end-5);
  underscore_pos = regexp(fiber,'_');
  if ~isempty(underscore_pos)
    fiber(underscore_pos-1:end) = [];  % handle files like  "20230821/eA_10°.txt"
  end
  
  COMfile = [file_full(1:slashes(end)),fiber,'COM.txt'];
  
  fid = fopen(COMfile);
  if fid == -1
    T = NaN;  % File not found
    return
  end
  c = textscan(fid,'%s','delimiter','\n');
  fclose(fid);
  lines = c{1};
  nlines = numel(lines);
  TB = [];
  for j = 1:nlines
    if contains(lines{j},'temperatureB')
      [~,pos] = regexp(lines{j},'temperatureB =');
      TB = [TB;str2double(lines{j}(pos+1:numel(lines{j})))];
    end
  end
  T = round(mean(TB,'omitnan'),2);
  % if isnan(T)
  %   T = 20;   % Just to choose a default
  % end
end