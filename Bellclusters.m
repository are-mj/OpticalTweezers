function [bellth,bellstd,n,rms] = Bellclusters(testcase,plotting)
% Exploring the effect of cluster definitions on Bell parameters
% and residuals (expressed at Root Mean Square residuals (rms))
% Figures are shown on extra screen.  Modify positions if needed

  if nargin < 2
    plotting = 0;
  end
  load filesets TPold TPk TPnew
  TPall = [TPold;TPk;TPnew];
  TPgood = remove_bad(TPall);


  % ********* Plot Basic shapes for building clusters:
  shapes = clustershapes(plotting);
  
  dF = 1;
  % Select table to analyse
  % tbl = TPall;
  tbl = TPgood;
  
  % ************* % Histograms of temperature and Fdot
  if plotting
    figure;
    histogram(tbl.Temperature,0.5:40.5)
    xlabel('°C')
    xline([10,20]);
    title('Temperature histogram')
    set(gcf,'position',[2492,216,560,420]);
    figure;
    histogram(log10(tbl.Fdot(tbl.Fdot>0)),100)
    set(gca,'XTick',[log10(0.5),0,log10(2),log10(5),1,log10(20),log10(50),2]);
    set(gca,'XTickLabel',{'0.5';'1';'2';'5';'10';'20';'50';'100'});
    xlim([log10(0.5),2]);
    xline([log10(2),log10(20)]);
    xlabel('Fdot (nm/s)')
    title('Histogram of force rate of change at rip')
    set(gcf,'position',[3052,216,560,420]);
  end
  
  % ************  Temperature grouping:
  lowT = tbl.Temperature<10;
  mediumT = tbl.Temperature >=10 & tbl.Temperature < 20;
  highT = tbl.Temperature>=20;

  % ************   Fdot grouping:
  slow = tbl.Fdot < 2;
  normal = tbl.Fdot >= 2 & tbl.Fdot < 20;
  fast = tbl.Fdot >= 20;

  % ****************** Select shapes for clusters ***************
  shapenos = {[2,3,5,6,8,9],4,1};
  switch testcase
    case 1
      % shapenos = {[2,3,5,6,8,9],[4],[1]};
      selected = lowT & normal;
      casetext = "T < 10°C";
    case 2
      % shapenos = {[2,3,5,6,8,9],[4]};
      % shapenos = {[1:9,11]};
      selected = mediumT & normal;
      casetext = "10°C <= T < 20°C";
    case 3
      % shapenos = {[1,2,3,4,5,6,8,9]};
      shapenos = {[1:9,11]};
      selected = highT & normal;
      casetext = "T >= 20°C"; 
  end
  title1 = strcat("Scatter plot with clusters. ",casetext);
  title2 = strcat("Best model fit. ",casetext);
  
  % *********** Select rips for clusters  ***************
  n_clusters = numel(shapenos);
  Cluster = clusters(tbl,shapenos,shapes,selected,plotting);
  if plotting
    set(gcf,'position',[1932,-282,560,420])
    title(title1)
  end
  
  Cluster = Cluster&selected;
  allclusters = any(Cluster,2);
  force = tbl.Force(allclusters&selected);
  [pdtot,Fbins] = probdens(force,dF);
  Ftot = (Fbins(1:end-1)+Fbins(2:end))/2;
  
  n = sum(Cluster);  % Number of rips in each cluster
  
  dx0 = 2;log10k0 = -3; theta0 = [dx0;log10k0];
  if plotting
  figure; hold on;
    for i = 1:n_clusters
      histogram(tbl.Force(Cluster(:,i)&selected),4:51);
    end
    box on
    set(gcf,'position',[2492,-282,560,420]);
  end
  Fplot = linspace(5,55);
  bellth = zeros(2,n_clusters);
  bellstd = zeros(2,n_clusters);
  pdbell = zeros(100,n_clusters);
  pdcalc = zeros(numel(Ftot),n_clusters);
  for i = 1:n_clusters
    [pd_obs,edges] = probdens(tbl.Force(Cluster(:,i)&selected),1);
    Tmean = mean(tbl.Temperature(Cluster(:,i)&selected));
    Fdotmean = mean(tbl.Fdot(Cluster(:,i)&selected));
    [bellth(:,i),bellstd(:,i)] = fit_Bell_unfold(pd_obs,edges,Tmean,Fdotmean,theta0);
    pdbell(:,i) = Bell_unfold_probability(bellth(:,i),Fplot,Tmean,Fdotmean);
    pdcalc(:,i) = Bell_unfold_probability(bellth(:,i),Ftot,Tmean,Fdotmean);
  end
  w = n/sum(n);
  residual = pdcalc*w'-pdtot;
  rms = sqrt(residual'*residual/numel(Ftot));
  if plotting
    figure; hold on;
    bar(Ftot,pdtot,dF)
    plot(Fplot,pdbell*w','r','LineWidth',1.5) 
    for i = 1:n_clusters
      plot(Fplot,pdbell(:,i)*w(i),'k');
    end
    set(gcf,'position',[3052,-282,560,420]);
    box on
    xlabel('Rip force (pN)');
    ylabel('Probability density (pN^-^1)');
    title(title2)
    legend('Observed','Calculated')
  end

  vars = ["x‡","log10(k0)"];
  fprintf('Case %d: %s\n',testcase,casetext);
  fprintf('Cluster:            %4d ',1);
  for j = 2:n_clusters
    fprintf('%11d ',j);
  end
  fprintf('\n');
  fprintf('Number of rips:     %4d',n(1));
  for j = 2:n_clusters
    fprintf('%11d ',n(j));
  end  
  fprintf('\n');
  for i = 1:2
    fprintf('%12s ',vars(mod(i+1,2)+1));
    for j = 1:numel(n)
      fprintf('%6.2f±%4.2f ', bellth(i,j),abs(bellstd(i,j)));
    end
    fprintf('\n')
  end
  fprintf('rms: %5.4f\n\n',rms);
end
