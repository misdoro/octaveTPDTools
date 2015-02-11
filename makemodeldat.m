#!/usr/bin/octave --persist
#Model a TPD dat file

par.bline=0;
par.scale=1;
par.ml=2e-8;
par=loadParamFile(par);
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
par.rate=useArgument(argv(),2,par.rate);
pkg load odepkg
pkg load general


tpd.mids=par.mid;
tpd.rate=par.rate;
maxt=(par.maxT-par.minT)/par.rate;
tpd.t=linspace(0,maxt,par.npT)';
tpd.T=linspace(par.minT,par.maxT,par.npT)';

#parv=[par.v,par.theta,par.monolay*par.rate,par.E];
#tpd.i=calcP1(tpd,parv)+par.bline;
#tpd.i=calcPn(tpd,parv,par)+par.bline;
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

##########################
#Output data truncation and save
if (nargin>=1);
	filename=useArgument(argv(),1,'output.dat');
	if (!isempty(strfind(filename,".dat" )))
		tpd.version=20140120;
		save("-binary",filename,"tpd","dose","press");
	endif
else
	printf("Usage: modelTPD filename rate");
endif
exit(0);