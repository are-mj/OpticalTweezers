function shortname = shorten_filename(fullname,levels)
% Utility that returns the last levels of a path and file combination
%   levels: the number of levels in shortname.  Default: 2
% Example: fullname = "C:\Users\are\Documents\20230712\aa.txt"
%          shorten_filename(fullname,2) returns "20230712/aA.txt"

  if nargin < 2
    levels = 2;
  end
  filename = string(strrep(fullname,'\','/'));  % Use Unix separator
  slashes = [0,regexp(filename,'\/')];  % Position of '/' in files
  fn = char(filename);  % Translate to character array
  maxlevels = numel(slashes);
  shortname = fn(slashes(end)+1:end);
  for i = 1:min(maxlevels,levels)-1
    shortname = [fn(slashes(end-i)+1:slashes(end-i+1)),shortname];
  end
  shortname = string(shortname);  % Convert back to string
