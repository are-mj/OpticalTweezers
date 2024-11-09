% model_fitting_deno.m
% Matlab script with the commands from: Fitting models to experiment.pptx

Dudko = false;  % True: Fit Bell and Dudko  models  fale: Fit Bell model only

% Replace the folder names with thise that are relevant in your case:
folders = ["02022022","02032022","02042022","02092022","02102022","02142022","02152022","02162022","02182022"];

files = textfile_names(folders);  % Return the names of all *.txt files in the foders, except *COM.txt files

% [TRIP,TZIP] = collect_tables(files);  % Comment out this is TRIP and TZIP have elready been created

% Select subset of results:
ok = TRIP.Temperature>20 & TRIP.Pullingspeed > 500;
force = TRIP.Force(ok);

dF = 1; % distance between force bin edges    
[pd_obs,edges,n_obs] = probdens(force,dF);
F_list = (edges(1:end-1)+edges(2:end))/2;     % force bin midpoints
figure; bar(F_list,pd_obs,dF);

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
title("Data for T>20\circC, Puling speed > 500 nm/s")
legend("Observed","Bell model","location","northwest")
if ~Dudko
  % Write results in figure
  text(8,0.12,sprintf('n_o_b_s = %d',n_obs))
  text(8,0.11,sprintf('x^‡ = %.2f ± %.2fnm',theta(1),theta_std(1)))
  text(8,0.10,sprintf('log10(k0) = %.2f ± %.2f',theta(2),theta_std(2)));
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
  text(8,0.11,sprintf('n_o_b_s = %d',n_obs))
  text(8,0.10,sprintf('ΔG^‡ = %.1f ± %.1fzJ',theta(1),theta_std(1)))
  text(8,0.09,sprintf('x^‡ = %.2f ± %.2fnm',theta(2),theta_std(2)))
  text(8,0.08,sprintf('log10(k0) = %.2f ± %.2f',theta(3),theta_std(3)))
end


