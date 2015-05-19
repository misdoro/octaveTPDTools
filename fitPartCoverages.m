function fitpar=fitPartCoverages(tpd,fitpar,\
fitopts=optimset("Display","final","MaxIter",3000,"TolX",1e-5),debug=0)
%Fit partial coverages of the energy distribution
%function fitpar=fitCoverages(tpd,fitpar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: For sub-ML coverages use scaling factor instead of resolving PW eqn every time,
%that should speed-up everything.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf("Optimising energy distribution\n");
fh=@(thetavar)fitcovs(tpd,fitpar,thetavar);

[thetao,ssq]=fminunc(fh,fitpar.thetas,fitopts);
fitpar.thetas=abs(thetao);
endfunction

function ssq=fitcovs(tpd,fitpar,thetas)
  fitpar.thetas=abs(thetas);
  dbg=fitpar.thetas
  p=modelTPDmc(tpd.T,fitpar);
	ssq=sumsq(p-tpd.i);
  ssqi=ssq;
  if (isfield(fitpar,'penalty'))
    ssq+=fitpar.penalty*penalty(fitpar.thetas,fitpar.dE);
  endif
  ssq*=1e20;
  if isfield(fitpar,'debug')
    printf("SSQ %e (+penalty=%e)\n",ssqi,ssq);
  endif
endfunction;

function result=penalty(thetas,dE,debug=0);
  result=sumsq(diff(thetas,2));
  #printf("penalty%e\n",result);
endfunction;
