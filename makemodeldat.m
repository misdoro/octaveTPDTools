#!/usr/bin/octave -q
#Model a TPD dat file

par.bline=0;
par.scale=1;
par.ml=2e-8;


#Load fit data if available
fit=loadFitFile();
if (isfield(fit,"fitdE"))
  fdE=fit.fitdE;
  par.v=fdE.v;
  par.E0=fdE.E0;
  par.dE=fdE.dE;
  par.thetas=fdE.thetas;
  par.rate=fdE.rate;
endif

#Load param file to override fit data
par=loadParamFile(par);

#Verify all required fields
fields={'mid','rate','npT','minT','maxT','v','E0','dE','thetas','bline'};
for i=1:length(fields)
  if (~isfield(par,fields(i)))
    printf("Param %s is missing",fields{i});
    printf("\n\
    param.m fields:\n\
    par.mid=mass\n\
    par.rate=K/sec\n\
    par.npT=number of temperature points\n\
    par.minT=min temperature\n\
    par.maxT=max temperature\n\
    par.v=prefactor\n\
    par.E0=initial energy\n\
    par.dE=energy step\n\
    par.thetas=[array of partial coverages]\n\
    par.bline=constant baseline\n\
    optional:\n\
    par.noise=1 to enable noise\n\
    par.snr=snr value\n\
    par.constnoise=constant noise\n\
    ");
    exit(0);
  endif
endfor

#Rate may be specified as second command-line parameter
par.rate=useArgument(argv(),2,par.rate*60)/60;
pkg load odepkg
pkg load general


tpd.mids=par.mid;
tpd.rate=par.rate;
maxt=(par.maxT-par.minT)/par.rate;
tpd.t=linspace(0,maxt,par.npT)';
tpd.T=linspace(par.minT,par.maxT,par.npT)';

tpd.i=modelTPDmc(tpd.T,par)+par.bline;
if(isfield(par,"noise")&& par.noise)
  snr=1000;
  if (isfield(par,"snr"))
    snr=par.snr;
  endif
  noiselvl=max(tpd.i)/snr
  if (isfield(par,"constnoise"))
    noiselvl+=par.constnoise;
  endif
  tpd.i=tpd.i.+stdnormal_rnd(par.npT,1).*noiselvl;
endif
tpd.iN=tpd.i;
press.t=tpd.t+0.1;
press.T=tpd.T;
press.To=press.T;
press.p=tpd.i*23;

tpd.pi=interp1(press.t,press.p,tpd.t);

figure(1);
hold on;
plot(press.T,press.p,"color","red",'linewidth',2);
plot(tpd.T,tpd.i,"color","green",'linewidth',2);

tpd.integral=trapz(tpd.t,tpd.i);
tpd.model=1;
tpd.par=par;

dose.integral=tpd.integral*0.1;
dose.T=par.minT;
dose.mids=tpd.mids;
dose.np=20;
dose.time=10;
dose.t=linspace(0,dose.time,dose.np)';
dose.i=dose.integral./10.*ones(1,dose.np)';
dose.iN=dose.i;


function printTPDInfo(par)
	printf("\nkey TPD parameters:\n\
          rate=%f\n\
          v=%e\n\
          E0=%f\n",par.rate,par.v,par.E0);
  if ((np=length(par.thetas))>1)
    printf("          Ei=[");
    printf(" %.3f ",linspace(par.E0,par.E0+(np-1)*par.dE,np));
    printf("]\n");
    
    printf("          thetas=[");
    printf(" %.3f ",par.thetas);
    printf("]\n");
    
    
  else
    printf("          theta=%f\n",par.thetas);
    
  endif
endfunction;


##########################
#Output data truncation and save
if (nargin>=1);
	filename=useArgument(argv(),1,'output.dat');
	if (!isempty(strfind(filename,".dat" )))
		tpd.version=20140120;
		save("-binary",filename,"tpd","dose","press");
		printf("Saved %s",filename);
		printTPDInfo(par);
	endif
else
	printf("Usage: modelTPD filename rate");
endif
exit(0);
