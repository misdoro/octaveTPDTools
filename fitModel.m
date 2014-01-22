function result=fitModel(mytpd,param,result,press,dose);
	result=plotTPD(mytpd,param,result);
	par.mids=mytpd.mids;
	
	
	par=loadParamFile(par);
	
	para=[par.v,dose.integral/par.monolayd,par.asc*par.monolay*mytpd.rate,par.E];
	ptotn=calcPn(mytpd,para,par);
	
	plot(mytpd.T,ptotn,'color','green');
	plot(mytpd.T,(mytpd.i-ptotn)*10,'color','black');
	
	
	plot(mytpd.T,mytpd.i,".");
	drawnow();
	input("hold on before fitting");
	
	
	fh=@(para)model(para,mytpd,par);
	options(1)=0;
	options(6)=2;
	[paro,mindiff]=fminsearch(fh,para,options,[]);
	paro
	
	popt=calcPn(mytpd,paro,par);
	plot(mytpd.T,popt,"color","red");
	plot(mytpd.T,(mytpd.i-popt)*10,'color','blue');
endfunction;


	
function ssq=model(parv,tpd,par)
	tic();
	p=calcPn(tpd,parv,par);
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