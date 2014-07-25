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
  
  fns=strsplit(mytpd.filename,".");
  fn=fns{1};
  if (isfield(par,fn))
    parfs=getfield(par,fn);
    cutT.min=parfs.minT;
    cutT.max=parfs.maxT;
    mytpd=cutTemp(mytpd,cutT);
  endif
    
  if (isfield(par,"decimate") && par.decimate)
    factor = round(length(mytpd.i)/par.np)
    mytpd=decimateTPD(mytpd,factor);
  endif
	dosecov=doseintg/par.monolayd;
  tpdcov=mytpd.intg/par.monolay
	para=[par.v,dosecov,par.monolay*mytpd.rate,par.E];
	ptotn=calcPn(mytpd,para,par);
	
  mytpd.i=mytpd.i-par.bline;
	plot(mytpd.T,ptotn,'color',mytpd.color);
	plot(mytpd.T,(mytpd.i-ptotn)*10,'color','black');
	ssq=fitcoverage((dosecov),para,mytpd,par);
  printf("Initial SSQ: %e\n",ssq);
	
	plot(mytpd.T,mytpd.i,".");
	drawnow();
  
	if (index(param.tools,'F'))
    for i=[1,2]
    printf("Coverage optimisation\n");
    fh=@(fp)fitcoverage(fp,para,mytpd,par);
    opts=optimset("Display","final","TolX",1e-9)
    fp=para(2);
  	[paro,mindiff]=fminsearch(fh,fp,opts);
    printf("Final SSQ: %e\n",mindiff);
    #paro= fminunc(fh,fp,opts);
    para(2)=paro(1);
  	popt=calcPn(mytpd,para,par);
  	#plot(mytpd.T,popt,"color","red");
	  #plot(mytpd.T,(mytpd.i-popt)*10,'color','blue');
    
    printf("Ea optimisation\n");
    fp=para(4);
    fh=@(fp)fitEa(fp,para,mytpd,par);
    [paro,mindiff]=fminsearch(fh,fp,opts);
    printf("Final SSQ: %e\n",mindiff);
    para(4)=paro(1);
  	popt=calcPn(mytpd,para,par);
    #plot(mytpd.T,popt,"color","green");
	  #plot(mytpd.T,(mytpd.i-popt)*10,'color','blue');
    endfor
    plot(mytpd.T,popt,"color","red");
	  plot(mytpd.T,(mytpd.i-popt)*10,'color','green');
    printf("Final values: Ea=%f eV, theta=%f ML\n",para(4),para(2));
  endif;
endfunction;

function ssq=fitcoverage(fp,parv,tpd,par)
  parv(2)=fp(1);
	p=calcPn(tpd,parv,par);
	ssq=sumsq(p-tpd.i);
endfunction;

function ssq=fitEa(fp,parv,tpd,par)
  parv(4)=fp(1);
	p=calcPn(tpd,parv,par);
	ssq=sumsq(p-tpd.i);
endfunction;