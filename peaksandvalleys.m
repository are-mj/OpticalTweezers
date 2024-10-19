function  [peakpos,valleypos] = peaksandvalleys(f,threshold,lim,plotting)
  % Find the  peaks and valleys of oscillating signal f
  %   The signal may be noisy and have varying oscillation frequency
  % Inputs:
  %          f:         Column vector with sinal series
  %          threshold: Separates f into high and low periods
  %          lim:       If max(f) < lim within a high f period, that perod
  %                     is classified as noise and marked as invalid.
  % Outputs:
  %          peaks:     Indices of maximum f within all valid high f periods 
  %          valleys:   Indices of minimum f within all low f periods

  if nargin < 4
    plotting = 0;
  end

  % Make sure no f values exactly match threshold:
  exact_hit = f==threshold;
  f(exact_hit) = f(exact_hit)*(1+5*eps);
  nf = size(f,1);

  % Use run length encoding to find lengths of all period of high or low f
  r = rle(f>threshold);
  if r(1,1)  % Starts with high value, so the first interval is low
    r = [[0;1],r];
  end
  ix = cumsum(r(2,:)); % index of start of high f periods
  nix = length(ix);
  high = [ix(1:2:nix-1)',ix(2:2:nix)'-1];
  fine = @(i,j) any(f(i:j)>lim);
  c = arrayfun(fine,high(:,1),high(:,2));
  high = high(c,:);

  pfun = @(i,j) peakfun(i,j,f);
  peakpos = arrayfun(pfun,high(:,1),high(:,2));
  vfun = @(i,j) valleyfun(i,j,f);
  % valleypos = arrayfun(vfun,[high(1:end-1,2)],[high(2:end,1)]);
  valleypos = arrayfun(vfun,[1;high(:,2)],[high(:,1);nf]);
  if plotting
    figure;
    plot(f);
    hold on;
    % plot(high(:,1),f(high(:,1)),'^r');
    % plot(high(:,2),f(high(:,2)),'vk');  
    plot(peakpos,f(peakpos),'ok');
    plot(valleypos,f(valleypos),'om')
  end
end


function peakpos = peakfun(i,j,f)
  [~,m] = max(f(i:j));
  peakpos = i+m-1;
end

function valleypos = valleyfun(i,j,f)
  [~,m] = min(f(i:j));
  valleypos = i+m-1;
end






