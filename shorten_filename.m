function shortname = shorten_filename(fullname,levels)
% Utility that returns the last levels of a path and file combination
% Example: fullname = "C:\Users\are\Documents\20230712\aa.txt"
% shorten_filename(fullname,2)
%   20230712/aA.txt
  filename = strrep(fullname,'\','/');  % Use Unix separator
  slashes = [0,regexp(filename,'\/')];  % Position of '/' in files
  fn = [char(filename)];  % Translate to character array
  maxlevels = numel(slashes);
  shortname = fn(slashes(end)+1:end);
  for i = 1:min(maxlevels,levels)-1
    shortname = [fn(slashes(end-i)+1:slashes(end-i+1)),shortname];
  end
  % if numel(slashes) <= levels
  %   shortname = fn;  % form: 'xx.txt' or '<subfolder>/xX.txt'
  % else
  %   shortname = fn(slashes(end-1)+1:end); % '<subfolder>/xX.txt'
  % end
  shortname = string(shortname);  % Convert back to string
