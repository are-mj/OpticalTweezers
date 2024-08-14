function [Tp,Tr,Tp_,Tr_] = compare(file)
% Compare original and new analyse_file
  home = [getenv('HOMEDRIVE') getenv('HOMEPATH'),'\OneDrive\Documents\'];
  newfiles = 'GitHub\OpticalTweezers';
  oldfiles = 'Nature_files';
  cd([home,newfiles]);
  par = par_single;
  [Tp,Tr] = analyse_experiment(file,1);
  ylim([0,60])
  ax1 = gca;
  set(gcf,"position",[1920,220,1901,420]);
  zoom on
  cd([home,oldfiles]);
  [Tp_,Tr_] = analyse_file(file,1);
  ylim([0,60])
  ax2 = gca;
  set(gcf,"position",[1920,-270,1901,420]);
  zoom on
  cd([home,newfiles]);
  linkaxes([ax1,ax2],"xy")
  fprintf('Unfolding: %d %d, Refolding: %d %d\n',height(Tp),height(Tp_),height(Tr),height(Tr_));
  pause;
end