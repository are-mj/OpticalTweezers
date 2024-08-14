function [section_borders,xx] = rle_smooth(x,brief)
% Smooths out brief periods of anomalours values in x
%
% Assumes x consists of perieonds repeated values, interspersed with
% brief periods of anomalous x values.  THis function changes the 
% anomolous entries to the surrounding values.
%  Inputs: x: 1D array of values
%  brief:  maximum length of anomalous periods
  
  x = x(:)';   % make sure x is row array
  r0 = rle(x);
  ix = cumsum(r0(2,:));
  short = r0(2,:)<brief;  % Identify anomalous periods
  r0(:,short)=[];
  r0(2,:) = diff([0,ix(~short)]);  % run lengths in reduced rle matrix
  r = rle(r0(1,:));    % runlemgths of values in reduces matrix
  ix(short)=[];        % Cumulative indices for reduced matrix only
  ixrr = cumsum(r(2,:));
  ixr = ix(ixrr);
  section_borders = [1,ixr];
  if nargout > 1
    r(2,:) = diff([0,ixr]);  
    % decode r:
    xx = [];
    for i = 1:numel(ixr)
      xx = [xx,ones(1,r(2,i))*r(1,i)];
    end
  end
end