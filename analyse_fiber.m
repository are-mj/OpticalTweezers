function [Trip,Tzip,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_fiber(fiber,plotting,par)
% Analyse all files for the same fiber in results table T
% fiber: e.g. '20230722/o' to match '20230722/oA', '20230722/oB' etc

if nargin < 3
  par = params;
end
if nargin < 2
  plotting = 0;
end
d = dir(fullfile(datafolder,[char(fiber),'*.txt']));
files = strings(numel(d),1);
for i = 1:numel(d)
  files(i) = fullfile(d(i).folder,string(d(i).name));
end
files(contains(files,"COM"))=[];

[Trip,Tzip,pull,relax,t,f,x,T,peakpos,valleypos] = analyse_experiment(files(1),0,par);
if numel(files) >1
  for i = 2:numel(files)
    [Trip0,Tzip0,pull0,relax0,t0,f0,x0,T0,peakpos0,valleypos0] = analyse_experiment(files(i),0,par);
    Trip = [Trip;Trip0];
    Tzip = [Tzip;Tzip0];
    pull = [pull;pull0];
    relax = [relax;relax0];
    peakpos = [peakpos;peakpos0+numel(t)];
    valleypos = [valleypos;valleypos0+numel(t)];  
    t = [t;t0];
    f = [f;f0];
    x = [x;x0];
    T = [T;T0];
  end
end
pt = vertcat(pull.time);
pf = vertcat(pull.force);
rt = vertcat(relax.time);
rf = vertcat(relax.force);

if plotting
  figure;
  plot(t,f,pt,pf,'*r',rt,rf,'ok')
end

