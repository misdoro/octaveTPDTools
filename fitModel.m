function result=fitModel(mytpd,param,result,press,dose);
	#result=plotTPD(mytpd,param,result,press,dose);
	par.mids=mytpd.mids;
	
	doseintg=0;
	if (isfield(mytpd,"version") && mytpd.version>=20140120 && ~mytpd.model)
		doseintg=calculateDoseIntegral(getMassData(dose,[],param.selectedmass));
	else
		doseintg=mytpd.doseintg;
	endif
	par=loadParamFile(par);
  if (isfield(par,"decimate"))
    factor = round(length(mytpd.i)/par.np)
    mytpd=decimateTPD(mytpd,factor);
  endif
	doseintg/par.monolayd
	para=[par.v,mytpd.intg/par.monolay,par.asc*par.monolay*mytpd.rate,par.E];
	ptotn=calcPn(mytpd,para,par);
	
  mytpd.i=mytpd.i-par.bline;
	plot(mytpd.T,ptotn,'color',mytpd.color);
	plot(mytpd.T,(mytpd.i-ptotn)*10,'color','black');
	
	
	plot(mytpd.T,mytpd.i,".");
	drawnow();
	if (index(param.tools,'F'))
    input("hold on before fitting");
  	fh=@(para)model(para,mytpd,par);
  	#options(1)=0;
  	#options(6)=2;
    #opts=optimset()
  	[paro,mindiff]=fminsearch(fh,para);
  	paro
	
  	popt=calcPn(mytpd,paro,par);
  	plot(mytpd.T,popt,"color","red");
	  plot(mytpd.T,(mytpd.i-popt)*10,'color','blue');
  endif;
endfunction;

function ssq=model(parv,tpd,par)
	tic();
	p=calcPn(tpd,parv,par);
	toc();
	ssq=sumsq(p-tpd.i)
endfunction;