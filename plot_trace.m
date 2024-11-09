function h = plot_trace(ax,s,r)
% Plot force vs. x trace curves with rips and zips
% Alternative calls:
%  h = plot_trace(ax,s,r);  % Will plot in axes ax
%  h = plot_trace(s,r);   % Will plot in current axes or create new figure
% Input:
%   ax:  Axis (e.g. app.UIAxes)  Default; gca (current axes)
%   s:   stretch trace struct
%   r:   relax trace struct (optional)
% Output: h: handle struct to stretch and relax trace plots
%   h(1): handle struct to stretch, h(2): handle struct to relax.
%     h(i).Graph: trace line
%     h(i).Rips: handles to rip/zip markers
%     h(i).Riplines: [n_rips×2 double] deltax, rip to previous deltax
%     h(i).Slope: [n_rips×2 double]  Slope lines before and after rip

% Author: Are Mjaavatten
% Version 2024-06-14
  
  % Handle the case whne no axes are supplied  (e.g. plot_trace(s,r))
  if isa(ax,'struct')
    if nargin > 1
      r = s;
    end    
    s = ax;
    ax = gca;
  end

  par = params;
  cla(ax);
  hold(ax,'on') 
  % Plot stretching trace:
  h(1) = ripplot(ax,s,1,par);

  if exist('r','var')
    % Plot relaxing trace:
    h(2) = ripplot(ax,r,-1,par);
  end
  xlabel(ax,'Trap position (nm)')
  ylabel(ax,'Force (pN)')
  title(ax,sprintf('%3.0fs ≤ time ≤ %3.0fs',s.t(1),s.t(end)),Interpreter="none")
  % Make sure that the force axis starts at 0:
  lims = ylim(ax);
  ylim(ax,[0,lims(2)])
  box(ax,'on');
  zoom(ax,'on');
end
  
function h = ripplot(ax,s,sgn,par)
% Plots trace graph for stretch or relax trace s.  Returns struct of graphics
% handles to the trace elements
%   sgn: 1: stretch, -1: relax

  % Make sure all fields are defined:
  h.Graph = [];
  h.Rips = [];
  h.Riplines = [];
  h.Slope = [];
  h.Ripno = [];
  h.totrip = [];

  graphcolor = 'r';
  ripcolor = 'b';
  if sgn < 0
    graphcolor = 'b';
    ripcolor = 'r';
  end
  h.Graph = plot(ax,s.x,s.f,graphcolor);
  if ~isfield(s,'rip_index')  % No rip or zip
    return
  end
  pos = s.rip_index;  % Short name to make code more readable
  n_rips = min(numel(pos),par.maxrips);
  for i = 1:n_rips
    % Rip marker:
    h.Rips(i) = plot(ax,s.ripx(i),s.force(i),['o',ripcolor],'MarkerSize',10);
    % deltax horizontal line:
    h.Riplines(i,1) = plot(ax,s.ripx(i)+[0,s.deltax(i)], ...
      s.force(i)*[1,1],'k'); 
    if i < n_rips 
      if sgn*(s.ripx(i)+s.deltax(i) - s.ripx(i+1)) > 0
        % if deltax line extends beyond next rip
        % draw extenion line from next rip to end of deltax line
        h.Riplines(i,2) = plot(ax,[s.ripx(i+1),s.ripx(i)+s.deltax(i)],...
          [s.force(i+1),s.force(i)],'k');
      end
    end
    % Plot linear fit before rip
    fitspan = round(numel(s.t)*par.maxfitfraction);
    brange = [max(1,pos(i)-fitspan),pos(i)];
    if i > 1  % Make sure brange does not extend beyond previous rip:
      brange = [max(brange(1),pos(i-1)),pos(i)-par.ripsteps];
    end
    h.Slope(i,1) = plot(ax,s.x(brange),polyval(s.pfx_b(i,:),s.x(brange)),'k','LineWidth',1);
    % Plot linear fit after rip.
    arange = [pos(i)+par.ripsteps,min(pos(i)+fitspan,numel(s.x))];
    if i < n_rips
      % to avoid confusing multiple lines, make sure arange does not 
      % overlap with brange for next rip
      arange = [arange(1),min(arange(end),pos(i+1)-par.ripsteps-par.linespan)];
    end    
    if diff(arange)>0 
      h.Slope(i,2) = plot(ax,s.x(arange),polyval(s.pfx_a(i,:),s.x(arange)),'k','LineWidth',1);
    end
    % h.Ripno(i) = text(ax,s.x(pos(i)),max(s.force(i)+sgn*3,0),num2str(sgn*i));
  end  
  % if sgn>0 & ~isempty(s.totrip)
  %   h.totrip = plot(ax,s.totrip.x,s.totrip.f,'-*k');
  % end

end