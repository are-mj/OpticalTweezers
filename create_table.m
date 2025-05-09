function T = create_table(st)
% Create Matlab results table array of pull or relax structs
% Input 
%    st: array of pull or relax structs
  if isempty(st)
    T = [];
    return
  end

  Filename = repmat(st(1).file,[numel(st),1]);
  Time = vertcat(st.time);
  Deltax = vertcat(st.deltax);
  Force = vertcat(st.force);
  Forceshift = vertcat(st.fstep);
  Trapx = vertcat(st.ripx);
  Fdot = vertcat(st.fdot);
  Slope = vertcat(st.slope);
  Pullingspeed = vertcat(st.pullingspeed); % Not implemented (yet)
  Temperature = vertcat(st.temperature);
  Timestep = vertcat(st.dt);
  Noise = vertcat(st.noise);
  % Lineno = vertcat(st.lineno);
  T = table(Filename,Time,Deltax,Force,Temperature,Forceshift,Trapx,...
          Fdot,Slope,Pullingspeed,Noise,Timestep);
end
