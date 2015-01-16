function fitpar=fitPartCoverages(tpd,fitpar,\
fitopts=optimset("Display","final","MaxIter",3000,"TolX",1e-5))
%Fit partial coverages of the energy distribution
%function fitpar=fitCoverages(tpd,fitpar)

fh=@(thetavar)fitcovs(tpd,fitpar,thetavar);

[thetao,ssq]=fminsearch(fh,fitpar.thetas,fitopts);
fitpar.thetas=thetao;
endfunction

function ssq=fitcovs(tpd,fitpar,thetas)
  fitpar.thetas=abs(thetas);
  p=modelTPDmc(tpd.T,fitpar);
	ssq=sumsq(p-tpd.i)
endfunction;