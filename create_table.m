function T = create_table(st)
% Create Matlab results table from pull or relax traces struct arrays
  if isempty(st)
    T = [];
    return
  end
  Filename = repmat(shorten_filename(st(1).file,2),[numel(st),1]);
  Time = vertcat(st.time);
  Deltax = vertcat(st.deltax);
  Force = vertcat(st.force);
  Forceshift = vertcat(st.fstep);
  Trapx = vertcat(st.ripx);
  Fdot = vertcat(st.fdot);
  Slope = vertcat(st.slope);
  Pullingspeed = vertcat(st.pullingspeed); % Not implemented (yet)
  Temperature = vertcat(st.temperature);
  % dt = vertcat(st.dt);
  % Lineno = vertcat(st.lineno);
  T = table(Filename,Time,Deltax,Force,Temperature,Forceshift,Trapx,...
          Fdot,Slope,Pullingspeed);
end