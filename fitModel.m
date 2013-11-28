function result=fitModel(mytpd,param,result);
	result=plotTPD(mytpd,param,result);
	mytpd.intg
	#par=[1.57e+14   1.2534e+00 2.02e-10  1.637e-01] #Xe on HOPG, low coverage
	#par=[5.55e+13   2.67e-10   1.396e+00   1.6403e-01 0.32 0.39]
	#par=[5.6209e+13 1.4259e+00   1.6403e-01   4.5996e-01   2.3794e-01]
	#par=[4.3017e+13   2.3249e+00   1.6461e-01   2.3508e-10]
	#par=[1e+14   5e+00   1.64e-01   1e-10]
	#par=[1e+13   4.7995e+00      7.3945e-11]
	#par=[3.7e+11   6.6528e+00   7.23e-11 0.136]#Xe multilayer on graphite
	#par=[3.7e+11   8.1e+00   7.23e-11   1.3600e-01]	#1.5A*s of Xe on graphite,
	#par=[3.7e+11   4.6e+00   7.23e-11   1.3600e-01]	#1A*s of Xe on graphite +60% coverage spread
	#par=[4.1418e+11   1.9488e+00   7.5530e-11   1.3819e-01]#Xe on graphite 5e-10 A*s
	#par=[4e+11   8   1.5904e-9   1.3578e-01]#Xe multi on ice
	#par=[ 5.21e+11   24   7.04e-11   1.36e-01]
	#par(1) is the v parameter (prefactor)
	#par(2) is a relative coverage
	#par(3) is a scaling multiplier
	#par(4) is the adsorption energy

	#par(2)=mytpd.intg/param.monolayer
	
	par=[8.8757e+11   5.4e+00   5.2895e-11   1.3781e-01]#Xe multi on graphite
	
	#[T,theta,p1]=modelTPD1(mytpd.T,par(2)+par(4),par(1),par(3),mytpd.rate);
	#[T,theta,p2]=modelTPD1(mytpd.T,par(2)-par(4),par(1),par(3),mytpd.rate);
	#pt1=p1*par(5)*2.67e-10;
	#pt2=p2*(1-par(5))*2.67e-10;
	
	#ptot=pt1+pt2;
	#plot(mytpd.T,pt1,'color','green');
	#plot(mytpd.T,pt2,'color','green');
	ptot=calcP1(mytpd,par);
	
	plot(mytpd.T,ptot,'color','red');
	
	tpdd=decimateTPD(mytpd,5);
	#tpdd=mytpd;
	
	ptotn=cumPress(tpdd,par);
	plot(tpdd.T,ptotn,'color','green');
	plot(tpdd.T,(tpdd.i-ptotn)*10,'color','black');
	
	
	plot(tpdd.T,tpdd.i,".");
	drawnow();
	
	#fh=@(par)model(par,tpdd);
	#options(1)=0;
	#options(6)=2;
	#[paro,mindiff]=fminsearch(fh,par,options,[]);
	#paro
	
	#popt=calcP1(mytpd,paro);
	#plot(mytpd.T,popt);
	
endfunction;
	
function press=calcP1(tpd,par)
	[T,theta,p0]=modelTPD1(tpd.T,par(2),par(1),par(4),tpd.rate);
	press=p0*par(3);
endfunction;

function ptotn=cumPress(tpd,par)
#Sum of N plots with varying coverage
	numsim=20;
	par(3)=par(3)/numsim;
	parm=repmat(par,numsim,1);
	cov=par(2);
	covspr=cov*0.4;
	#covspr=min(cov*0.6,2);
	printf("coverage spread: %3.2f\n",covspr);
	parm(:,2)=linspace(par(2)-covspr,par(2)+covspr,numsim);



	ptotn=zeros(length(tpd.T),1);
	for i=1:numsim
		ptotn+=calcP1(tpd,parm(i,:));
	endfor
endfunction;
	
function ssq=model(par,tpd)
	par
	tic();
	p=calcP1(tpd,par);
	toc();
	ssq=sumsq(p-tpd.i)
endfunction;
	
function ssq=model2(par,tpd)
	par
	tic();
	[T,theta,p1]=modelTPD1(tpd.T,par(2)+par(4),par(1),par(3),tpd.rate);
	[T,theta,p2]=modelTPD1(tpd.T,par(2)-par(4),par(1),par(3),tpd.rate);
	ptot=(p1*par(5)+p2*(1-par(5)))*2.67e-10;
	toc();
	p=ptot;
	ssq=sumsq(p-tpd.i)
endfunction;