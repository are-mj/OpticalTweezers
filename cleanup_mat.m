function cleanup_mat(matfile)
  if exist(matfile,"file")
    load(matfile,'Trip','Tzip');
    if any("Fitrange" == string(Trip.Properties.VariableNames))
      return
    end
    ripnos = height(Trip);
    if ripnos > 0
      Trip.Filename = repmat(shorten_filename(Trip.Filename(1)),[ripnos,1]);
      newripcols = zeros(ripnos,2);
      Trip = addvars(Trip,newripcols,'After',"Noise",'NewVariableNames',"Fitrange");
    end
    zipnos = height(Tzip);
    Tzip.Filename = repmat(shorten_filename(Tzip.Filename(1)),[zipnos,1]);
    if zipnos > 0
      newzipcols = zeros(height(Tzip));
      Tzip = addvars(Tzip,newzipcols,'After',"Noise",'NewVariableNames',"Fitrange");    
    end
    save(matfile,"Trip","Tzip");
  end