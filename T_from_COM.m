function T = T_from_COM(file)
% Read temperature from the COM file corresponding to the experiment file
%   Backup solution if the temperature is not coded in the staus variable
%   2023-10-16: Modified the construction of COM file name.  Now works for
%   abA.txt (but no longer for aAB.txt)
%   2023-12-07: Read only TemperatureB (as recommended by Steve B. Smith)
%   2024-10-29: Made extracting fiber name from file name more robust
%               'file' must include full path
  
  % Extract fiber name from file name
  [path,name,ext] = fileparts(file);
  name = char(name);  % To allow accessing individual characters 
  corename = name(isstrprop(name,'alpha'));
  if length(corename) >= 2 && isequal(isstrprop(corename(1:2),'lower'),[1 0]) % aBxxx
    fiber = corename(1);
  elseif length(corename) >= 3 && isequal(isstrprop(corename(1:3),'lower'),[1 1 0]) % abCxxx
    fiber = corename(1:2);
  else
    error('Unable to extract fiber name from %s',name);
  end
  COMfile = fullfile(path,strcat(fiber,'COM.txt'));
  
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