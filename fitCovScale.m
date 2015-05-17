function fitpar=fitCovScale(tpd,fitpar,
fitopts=optimset("Display","final","MaxIter",3000,"TolX",1e-5))
%Fit coverage scale parameter
printf("Optimising coverage scaling\n");
fh=@(covsc)fitcovsc(tpd,fitpar,covsc);

[covo,ssq,info]=fminunc(fh,fitpar.scale,fitopts);
info
fitpar.scale=abs(covo);
endfunction

function ssq=fitcovsc(tpd,fitpar,covsc)
  fitpar.scale=abs(covsc);
  p=modelTPDmc(tpd.T,fitpar);
	ssq=sumsq(p-tpd.i)*1e20;
  if isfield(fitpar,'debug')
    printf("SSQ %e\n",ssq);
  endif
endfunction;