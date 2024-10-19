function [Tp,Tr,Tp_,Tr_] = compare(file)
% Compare original and new analyse_file
  home = [getenv('HOMEDRIVE') getenv('HOMEPATH'),'\OneDrive\Documents\'];
  newfiles = 'GitHub\OpticalTweezers';
  oldfiles = 'Nature_files';

  cd([home,oldfiles]);
  [Tp_,Tr_] = analyse_file(file,1);
  ylim([0,60])
  % ax1 = gca;
  set(gcf,"position",[1920,220,1901,420]);
  zoom on  
  hold on;

  cd([home,newfiles]);
  % par = parameter_struct;
  [Tp,Tr] = analyse_experiment(file);
  plot(Tp.Time,Tp.Force,'or','MarkerSize',12);
  if ~isempty(Tr)
    plot(Tr.Time,Tr.Force,'.r');
  end
  % ylim([0,60])
  % ax2 = gca;
  % set(gcf,"position",[1920,-270,1901,420]);
  % zoom on
  % linkaxes([ax1,ax2],"xy")
  fprintf('Unfolding: %d %d, Refolding: %d %d\n',height(Tp),height(Tp_),height(Tr),height(Tr_));
  pause;
end