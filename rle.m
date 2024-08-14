function r = rle(x)
% Run length encoding of sequence x
% Adapted from Mohsen Nosratinia's answer at
%  https://stackoverflow.com/questions/12059744/run-length-encoding-in-matlab
  x = x(:)';  % make sure x is a row vector
  J=find(diff([x(1)-1, x]));  % Index where x changes value
  r=[x(J); diff([J, numel(x)+1])];
end