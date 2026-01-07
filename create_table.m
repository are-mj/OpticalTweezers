function T = create_table(st)
% Create Matlab results table array of pull or relax structs
% Input 
%    st: pull or relax struct

  if isempty(st)
    T = [];
    return
  end
  T = [];
  for k = 1:length(st)
    nrips = length(st(k).force);
    for i = 1:nrips
      % Filename = st(k).file;
      Filename = erase(st(k).file,datafolder+"/");  % remove datafolder part
      Time = st(k).time(i);
      Deltax = st(k).deltax(i);
      Force = st(k).force(i);
      Forceshift = st(k).fstep(i);
      Trapx = st(k).ripx(i);
      Fdot = st(k).fdot(i);
      Slope_b = st(k).pfx_b(i,1);
      Slope_a = st(k).pfx_a(i,1);
      Pullingspeed = st(k).pullingspeed(i); 
      Temperature = st(k).temperature(i);
      Topforce = st(k).topforce(i);
      Timestep = st(k).dt(i);
      Noise = st(k).noise(i);
      Fitrange = st(k).fitrange(i,:);
      Cycleno = st(k).cycleno(i);
      Work = st(k).work(i);
      T = [T;table(Filename,Time,Deltax,Force,Temperature,Forceshift,Trapx,...
              Fdot,Slope_b,Slope_a,Pullingspeed,Topforce,Noise,Fitrange,...
              Cycleno,Work,Timestep)];
    end
  end
end
