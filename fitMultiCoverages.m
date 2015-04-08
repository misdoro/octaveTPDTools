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
      fitpar.rate=mytpd.rate
    else
      fitpar=myInitFitPar(mytpd,param)
    endif
    if (isfield(param,"fitpenalty"))
      fitpar.penalty=param.fitpenalty;
    endif
    if (isfield(param,"fitmaxiter"))
      maxiter=param.fitmaxiter;
    else
      maxiter=1500;
    endif
    if (isfield(param,"debug"))
      fod="iter"
    else
      fod="final"
    endif
    fitopts=optimset("Display",fod,"MaxIter",maxiter,"TolX",1e-5);
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


function fitpar=myInitFitPar(mytpd,param)
  fit=loadFitFile();
  fitpar=struct();
  if (isfield(fit,"fitdE"))
    printf("Using best fit prefactor from fit.par\n");
    %fitpar=initFitParam(mytpd,param,)
    fitpar=fit.fitdE;
    fitpar.E0-=2*fitpar.dE;
    fitpar.dE=estimdE(mytpd,fitpar);
    fitpar.np=estimNP(mytpd,fitpar);
    fitpar.thetas=0.01*ones(fitpar.np,1);
  else
    fitpar=initFitParam(mytpd,param)
    fitpar.E0-=6*fitpar.dE;
    fitpar.thetas=0.01*ones(fitpar.np+6,1);
  endif
endfunction
    