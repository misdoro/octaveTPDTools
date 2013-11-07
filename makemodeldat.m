#!/usr/bin/octave --persist
#Model a TPD dat file

par.mids=[0 130]

par.rate=3/60;
par.np=300;
par.minT=30;
par.maxT=80;
par.theta=0.5;
par.v=1e14;
par.E=0.160;
par.monolay=2e-9;

pkg load odepkg

tpd.mids=par.mids;
tpd.rate=par.rate;
tpd.t=linspace(0,par.np*par.rate,par.np)';
tpd.T=linspace(par.minT,par.maxT,par.np)';


[To,th,p]=modelTPD1(tpd.T,par.theta,par.v,par.E,tpd.rate);
tpd.i=p*par.monolay*par.rate;

tpd.iN=tpd.i;
press.t=tpd.t+0.1;
press.T=tpd.T;
press.p=tpd.i*23;

tpd.pi=interp1(press.t,press.p,tpd.t);

figure(1);
hold on;
plot(press.T,press.p,"color","red",'linewidth',2);
plot(tpd.T,tpd.i,"color","green",'linewidth',2);

tpd.integral=trapz(tpd.t,tpd.i);
tpd.model=1;

dose.integral=tpd.integral*0.1;
dose.T=par.minT;

##########################
#Output data truncation and save
if (nargin>=1);
	filename=argv(){1};
	if (!isempty(strfind(filename,".dat" )))
		tpd.version=20131025;
		save("-binary",filename,"tpd","dose","press");
	endif
endif

