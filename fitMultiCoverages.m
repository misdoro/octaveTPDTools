function result=fitMultiCoverages(mytpd,param,result);
  [f.info, f.err, f.msg]=stat("fitcov.par");
  if (f.err>=0);
    load "fitcov.par";
    printf("Plotting saved fit data\n");
    fitcov=fitcov{mytpd.idx};
  else
    if (isfield(result,"fitcov"))
      printf("Using parameters from previous fit iteration\n");
      fitpar=result.fitcov{mytpd.idx-1};
    else
      fit=loadFitFile();
      fitpar=struct();
      if (isfield(fit,"fitdE"))
        printf("Using parameters from fit.par\n");
        fdE=fit.fitdE;
        fitpar=fit.fitdE;
        fitpar.E0-=2*fitpar.dE;
        fitpar.dE/=2;
        fitpar.thetas=0.01*ones(35,1);
      else
        fitpar.v=param.v;
        fitpar.E0=param.E0;
        fitpar.dE=param.dE;
        fitpar.thetas=param.thetas;
      endif
    endif
    if (isfield(param,"fitpenalty"))
      fitpar.penalty=param.fitpenalty;
    endif
    if (isfield(param,"fitmaxiter"))
      maxiter=param.fitmaxiter;
    else
      maxiter=1500;
    endif
    fitpar.rate=mytpd.rate
    fitopts=optimset("Display","iter","MaxIter",maxiter,"TolX",1e-5);
    fitcov=fitPartCoverages(mytpd,fitpar,fitopts);
  endif
  np=length(fitcov.thetas);
  Epts=linspace(fitcov.E0,fitcov.E0+(np-1)*fitcov.dE,np);
 
  figure(getFigIndex("covsites"));
  plot(Epts,fitcov.thetas,"linewidth",2,"color",mytpd.color);
  
  p=modelTPDmc(mytpd.T,fitcov);
  figure(getFigIndex("covfits"));
  plot(mytpd.T,mytpd.i,"color",mytpd.color);
  plot(mytpd.T,p,"color",mytpd.color,"linestyle","--");
  drawnow();
  #input("proceed to next figure");
  result.fitcov{mytpd.idx}=fitcov;
endfunction;