%Estimate dE of the fit
function fitdE=estimdE(mytpd,fit)
  %dummy fit parameters 
  fit.dE=0.1;
  fit.thetas=0.1;
  %Define our own dense temperature grid
  Ts=linspace(min(mytpd.T),max(mytpd.T),200);
  moddat=modelTPDmc(Ts,fit);
  #Find dE such that Tmax(tpd(E-dE)) is at 0.7*max(tpd(E)) on the onset
  iend=min(find(moddat>=(max(moddat)*0.7)));
  fitdE=abs(fit.E0-estimEaVTm(fit.v,Ts(iend),fit.rate));
endfunction