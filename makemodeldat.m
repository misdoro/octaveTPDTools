#!/usr/bin/octave --persist
#Model a TPD dat file

par.mids=[84];

par.rate=15/60;
par.np=321;
par.minT=25;
par.maxT=55;
par.theta=0.16;
par.thetas=0.001;
par.v=3.7e+11;
par.E=0.11;
par.Es=0.007;
par.monolay=2e-8;
par.bline=1e-13;
par.numsim=10;
par.ls=1;
par.parallel=1;


par=loadParamFile(par);

pkg load odepkg
pkg load general


tpd.mids=par.mids;
tpd.rate=par.rate;
maxt=(par.maxT-par.minT)/par.rate;
tpd.t=linspace(0,maxt,par.np)';
tpd.T=linspace(par.minT,par.maxT,par.np)';

parv=[par.v,par.theta,par.monolay*par.rate,par.E];
#tpd.i=calcP1(tpd,parv)+par.bline;
tpd.i=calcPn(tpd,parv,par)+par.bline;
if(isfield(par,"noise")&& par.noise)
  snr=1000;
  if (isfield(par,"snr"))
    snr=par.snr;
  endif
  maxisnr=max(tpd.i)/snr
  tpd.i=tpd.i.+stdnormal_rnd(par.np,1).*maxisnr;
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
	filename=argv(){1};
	if (!isempty(strfind(filename,".dat" )))
		tpd.version=20140120;
		save("-binary",filename,"tpd","dose","press");
	endif
else
	printf("Usage: modelTPD filename")
endif

exit(0);