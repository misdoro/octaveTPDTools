function fit=initFitParam(mytpd,param,v=param.v)
  fit.debug=param.debug;
  fit.ml=param.monolayer;
  fit.v=v;
  fit.rate=mytpd.rate;
  
  fit.scale=1;
  fit.E0=estimE0(mytpd,fit);
  fit.dE=estimdE(mytpd,fit);
  fit.np=estimNP(mytpd,fit);
  fit.thetas=0.1;
endfunction