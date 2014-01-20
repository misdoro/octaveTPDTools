#!/usr/bin/octave --persist
#Model a TPD dat file

par.mids=[0 84];

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
par.parallel=3;


par=loadParamFile(par);

pkg load odepkg
pkg load general

function press=calcP1o(tpd,parv)
	[T,theta,p0]=modelTPD1(tpd.T,parv(2),parv(1),parv(4),tpd.rate);
	press=p0*parv(3);
endfunction;

function press=calcP1ao(para)
	[T,theta,p0]=modelTPD1(para.T,para.b,para.a,para.d,para.rate);
	press=p0*para.c;
endfunction;

function press=calcP1(tpd,parv)
	[T,theta,p0]=modelTPDlsode(tpd.T,parv(2),parv(1),parv(4),tpd.rate);
	press=p0*parv(3);
endfunction;

function press=calcP1a(para)
	[T,theta,p0]=modelTPDlsode(para.T,para.b,para.a,para.d,para.rate);
	press=p0*para.c;
endfunction;

function ptotn=calcPn(tpd,parv,par)
	#Sum of tpds with variable parameters
	numsim=par.numsim;
	pars=parv;
	pars(3)=parv(3)/numsim;
	parm=repmat(pars,numsim,1);
	espr=par.Es;
	cspr=par.thetas;
	parm(:,4)=linspace(pars(4)-espr,pars(4)+espr,numsim);
	parm(:,2)=linspace(pars(2)-cspr,pars(2)+cspr,numsim);
	
	
	if (par.parallel<2);
		tic();
		ptotn=zeros(length(tpd.T),1);
	
		for i=1:numsim
			if (par.ls==1);
				ptotn+=calcP1(tpd,parm(i,:));
			else
				ptotn+=calcP1o(tpd,parm(i,:));
			endif
		endfor
		toc()
	else
	
		tic();
		for i=1:numsim
			para(i).a=parm(i,1);
			para(i).b=parm(i,2);
			para(i).c=parm(i,3);
			para(i).d=parm(i,4);
			para(i).T=tpd.T;
			para(i).rate=tpd.rate;
		endfor
		if (par.ls==1);
			ptotn=sum(pararrayfun(par.parallel,@calcP1a,para),2);
		else
			ptotn=sum(pararrayfun(par.parallel,@calcP1ao,para),2);
		endif;
		toc()
	endif;
endfunction;

tpd.mids=par.mids;
tpd.rate=par.rate;
maxt=(par.maxT-par.minT)/par.rate;
tpd.t=linspace(0,maxt,par.np)';
tpd.T=linspace(par.minT,par.maxT,par.np)';

parv=[par.v,par.theta,par.monolay*par.rate,par.E];
#tpd.i=calcP1(tpd,parv)+par.bline;
tpd.i=calcPn(tpd,parv,par)+par.bline;

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