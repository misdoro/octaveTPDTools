function fitpar=fitE0(tpd,fitpar,
fitopts=optimset("Display","final","MaxIter",3000,"TolX",1e-5))
%Fit initial energy point
printf("Optimising initial energy point\n");
fh=@(E0)fitcovsc(tpd,fitpar,E0);

[E0o,ssq]=fminsearch(fh,fitpar.E0,fitopts);

printf("Optimized value:%f\n",abs(E0o));
fitpar.E0=abs(E0o);
endfunction

function ssq=fitcovsc(tpd,fitpar,E0)
  fitpar.E0=abs(E0);
  p=modelTPDmc(tpd.T,fitpar);
	ssq=sumsq(p-tpd.i);
  if isfield(fitpar,'debug')
    printf("SSQ %e\n",ssq);
  endif
endfunction;