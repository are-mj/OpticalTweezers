function [T,instrument] = T_from_COM(file)
% Read temperature from the COM file corresponding to the experiment file
%   2023-10-16: Modified the construction of COM file name.  Now works for
%   abA.txt (but no longer for aAB.txt)
%   2023-12-07: Read only TemperatureB (as recommended by Steve B. Smith)
%   2024-10-29: Made extracting fiber name from file name more robust
%               'file' must include full path
  
  % Extract fiber name from file name
  [path,name] = fileparts(file);
  name = char(name);  % To allow accessing individual characters 
  corename = name(isstrprop(name,'alpha'));
  if length(corename) >= 3 && isequal(isstrprop(corename(1:3),'lower'),[1 1 0]) % abCxxx
    fiber = corename(1:2);
  else
    fiber = corename(1);
  end
  COMfile = fullfile(path,strcat(fiber,'COM.txt'));
  
  fid = fopen(COMfile);
  if fid == -1
    T = NaN;  % File not found
    warning('%s not found. Cannot determine temperaturw',COMfile);    
    return
  end
  c = textscan(fid,'%s','delimiter','\n');
  fclose(fid);
  lines = c{1};
  % Find lines containing 'temperatureB' or 'Instrument name'
  nlines = numel(lines);
  TB = [];
  instrument = '';
  for j = 1:nlines
    if contains(lines{j},'temperatureB')
      [~,pos] = regexp(lines{j},'temperatureB =');
      TB = [TB;str2double(lines{j}(pos+1:numel(lines{j})))];
    end
    if contains(lines{j},'Instrument name')
      [~,pos] = regexp(lines{j},'Instrument name =');
      instrument = strtrim(lines{j}(pos+1:numel(lines{j})));
    end    
  end
  T = round(mean(TB,'omitnan'),2);
end