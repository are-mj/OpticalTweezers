function [peakpos,valleypos,varargout] = decim(peakpos,valleypos,par,varargin)
% Reduce excessive number of data points in a set of time series.
%  This often makes it easier to detect rips nd zips.
%  Time series are smoothed and downsampled so that the number of data
%  points per trace is between 0.5 and 1 times par.maxpointspertrace.
% Used by analyse_experiment.m
% Typical use:
%   [peakpos,valleypos,t,f,x,T] = decim(peakpos,valleypos,t,f,x,T,par)
  if par.maxpointspertrace == 0
    varargout = varargin;
    return
  end
  idealpointspertrace = par.maxpointspertrace/2;
  pointspertrace = mean(diff(peakpos))/2;
  nvars = length(varargin);
  varargout = cell(1,nvars);
  if pointspertrace > 2*idealpointspertrace
    factor = round(pointspertrace/idealpointspertrace);
    nvars = length(varargin);  
    for i = 1:nvars
      smoothedvar = movmean(varargin{i},factor);
      varargout{i} = downsample(smoothedvar,factor);
    end
    peakpos = round(peakpos/factor);
    valleypos = round(valleypos/factor);
  else
    varargout = varargin;
  end
end