function ret=estimE0(mytpd,fit)
  [maxi,maxii]=max(mytpd.i);
  maxT=mytpd.T(maxii);
  ret=estimEaVTm(fit.v,maxT,mytpd.rate);
endfunction