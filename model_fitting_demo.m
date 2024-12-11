% model_fitting_deno.m
% Matlab script with the commands from: Fitting models to experiment.pptx

Dudko = false;  % True: Fit Bell and Dudko  models  fale: Fit Bell model only

% copy the contents of Summer school 2024/Testdata.zip to your Datafolder
folders = ["20230721","20230722","20230724","20230725"];

files = textfile_names(folders);  % Return the names of all *.txt files in the foders, except *COM.txt files
[TRIP,TZIP] = collect_tables(files);  % Comment out this if TRIP and TZIP have elready been created

% Select subset of results:
% ok = TRIP.Temperature<5 & TRIP.Pullingspeed < 200;
ok = TRIP.Temperature>5 & TRIP.Fdot < 15;

force = TRIP.Force(ok);

dF = 2; % distance between force bin edges    
[pd_obs,edges,n_obs] = probdens(force,dF);
F_list = (edges(1:end-1)+edges(2:end))/2;     % force bin midpoints
figure; bar(F_list,pd_obs,1);

Tmean = mean(TRIP.Temperature(ok),"omitmissing");
Fdotmean = mean(TRIP.Fdot(ok),"omitmissing");

% Fit Bell model:
dx0 = 2;lgk0 = -3; theta0 = [dx0;lgk0];
[theta,theta_std] = fit_Bell_unfold(pd_obs,edges,Tmean,Fdotmean,theta0);
fprintf("x‡ = %.2fnm ± %.2g\n",theta(1),theta_std(1))
fprintf("lg(k0) = %.2f  ± %.2g\n",theta(2),theta_std(2))

F_plot = linspace(edges(1),edges(end));
pd = Bell_unfold_probability(theta,F_plot,Tmean,Fdotmean);
hold on;
plot(F_plot,pd,'r','linewidth',2);
xlabel("Force (pN)")
ylabel("Probability density (pN^-^1)")
title("Data for T<20\circC, Puling speed < 200nm/s")
legend("Observed","Bell model","location","northeast")
limx = xlim;
limx(1) = min(limx(1),10);
xlim(limx);   % Make sure the x axis includes 10
limy = ylim;
linespace = (limy(2)-limy(1))/15;
if ~Dudko
  % Write results in figure
  text(10.5,limy(2)-linespace*1,sprintf('n_o_b_s = %d',n_obs))
  text(10.5,limy(2)-linespace*2,sprintf('x^‡ = %.2f ± %.2fnm',theta(1),theta_std(1)))
  text(10.5,limy(2)-linespace*3,sprintf('log10(k0) = %.2f ± %.2f',theta(2),theta_std(2)));
end

if Dudko
  % Fit Dudko model:
  dG0 = 80;dx0 = 2;lgk0 = -3; theta0 = [dG0;dx0;lgk0];
  par.nu = 1/2;par.model="DHS";
  [theta,theta_std,resnorm] = fit_Dudko_unfold(pd_obs,edges,Tmean,Fdotmean,theta0,par);
  fprintf("G‡ = %.1f±%.1f\n",theta(1),theta_std(1))
  fprintf("x‡ = %.2fnm ± %.2g\n",theta(2),theta_std(2))
  fprintf("lg(k0) = %.2f  ± %.2g\n",theta(3),theta_std(3))
  F_plot = linspace(edges(1),edges(end));

  hold on;
  plot(F_plot,pd,'--b','linewidth',2);
  xlabel("Force (pN)")
  ylabel("Probability density (pN^-^1)")
  title("Data for T>20\circC, Puling speed > 500 nm/s")
  legend("Observed","Bell model","Dudko model","location","northwest")

  % Write resukts in figure
  text(10.5,limy(2)-linespace*1,sprintf('n_o_b_s = %d',n_obs))
  text(10.5,limy(2)-linespace*2,sprintf('ΔG^‡ = %.1f ± %.1fzJ',theta(1),theta_std(1)))
  text(10.5,limy(2)-linespace*3,sprintf('x^‡ = %.2f ± %.2fnm',theta(2),theta_std(2)))
  text(10.5,limy(2)-linespace*4,sprintf('log10(k0) = %.2f ± %.2f',theta(3),theta_std(3)))
end


