function mark_cycle(ax,k,t,peakpos,valleypos,pullfirst)
% Mark cycle no k

  if isfield(ax.UserData,'region')
    delete(ax.UserData.region)  % remove existing region
  end

  peakfirst = valleypos(1)>peakpos(1); 
  if pullfirst
    if k < numel(valleypos)
      t1 = t(valleypos(k));
      t2 = t(valleypos(k+1));
    else
      return;
    end
  else
    if k < numel(peakpos)
      t1 = t(peakpos(k));
      t2 = t(peakpos(k+1));
    else
      return;
    end
  end
  region = xregion(ax,t1,t2);
  ax.UserData.region = region;
end
