function result=fitModel(mytpd,param,result);
	result=plotTPD(mytpd,param,result);
	mytpd.intg
	#v,mult,theta,Eads
	#par=[9.487e13, 1.613e-10, 2.06, 0.164]
	par=[1.1215e+14   1.1996e-10   1.2534e+00   1.6469e-01]
	#par=[9.9578e+13   2.1819e-11   1.3563e+00   1.6497e-01]
	#par=[1.08e14, 1.28e-10, 0.5, 0.163]#Xe
	#par=[5e14,1e-9,1,0.09]#Ar
	[T,theta,p]=modelTPD1(mytpd.T,par(3),par(1),par(4),mytpd.rate);
	plot(mytpd.T,p*par(2),'color','red');
	drawnow();
	
	tpdd=decimateTPD(mytpd,10);
	fh=@(par)model(par,tpdd);
	options(1)=1;
	options(6)=2;
	[paro,mindiff]=fminsearch(fh,par,options,[]);
	paro
	
	[T,theta,p]=modelTPD1(mytpd.T,paro(3),paro(1),paro(4),mytpd.rate);
	plot(mytpd.T,p*paro(2));
	
endfunction;
	
	
function ssq=model(par,tpd)
	par
	tic();
	[T,theta,p]=modelTPD1(tpd.T,par(3),par(1),par(4),tpd.rate);
	toc()
	p=p*par(2);
	ssq=sumsq(p-tpd.i)
endfunction;