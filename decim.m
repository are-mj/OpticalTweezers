function [data,factor] = decim(indata,peakpos,par)
% Reduce excessive number of data points in a set of time series.
%  This often makes it easier to detect rips nd zips.
%  Time series are smoothed and downsampled so that the number of data
%  points per trace is between 0.5 and 1 times par.maxpointspertrace.
% Used by analyse_experiment.m
% Typical use:
%   [data,factor] = decim(data,peakpos,valleypos,par)
  data = indata;
  factor = 1;
  idealpointspertrace = par.maxpointspertrace/2;
  pointspertrace = mean(diff(peakpos))/2;
  if pointspertrace > 2*idealpointspertrace
    factor = round(pointspertrace/idealpointspertrace); 
    data = movmean(indata,factor);
    data = data(round(factor/2):factor:end,:);
  end
end