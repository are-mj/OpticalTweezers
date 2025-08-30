function T = create_table(st)
% Create Matlab results table array of pull or relax structs
% Input 
%    st: pull or relax struct array
  if isempty(st)
    T = [];
    return
  end

  ncycles = length(st);  % Number of cycles with rips/zips
  T = [];
  for i = 1:ncycles
    nrips = length(st(i).force);  % Rips in cycle i
    for j = 1:nrips
      Filename = string(st(i).file);
      Time = st(i).time(j);
      Deltax = st(i).deltax(j);
      Force = st(i).force(j);
      Forceshift = st(i).fstep(j);
      Trapx = st(i).ripx(j);
      Fdot = st(i).fdot(j);
      Slope_b = st(i).pfx_b(j,1);
      Slope_a = st(i).pfx_a(j,1);
      Pullingspeed = st(i).pullingspeed(j); 
      Temperature = st(i).temperature(j);
      Topforce = st(i).topforce(j);
      Timestep = st(i).dt(j);
      Noise = st(i).noise(j);
      Fitrange = st(i).fitrange;
      Cycleno = st(i).cycleno(j);
      Work = st(i).work(j);
      T = [T;table(Filename,Time,Deltax,Force,Temperature,Forceshift, ...
        Trapx,Fdot,Slope_b,Slope_a,Pullingspeed,Topforce,Noise, ...
        Fitrange,Cycleno,Work,Timestep)];
      end
  end
end
