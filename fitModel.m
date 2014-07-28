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
    if isfield(par,"vs")
      vs=par.vs
    else
      vs=par.v
    endif
    fitret=[];
    for v=vs
      para(1)=v
      printf("Iteration for v=%e\n",v);
      ssq=0;
      for i=[1,2]
        printf("Coverage optimisation\n");
        
        #opts=optimset("Display","iter","TolX",1e-9);
        opts=optimset("Display","final","TolX",1e-9);
        
        [para,ssq]=optimEa(opts,para,mytpd,par);
        #plotoptimres(mytpd,para,par);
        
        [para,ssq]=optimCoverage(opts,para,mytpd,par);
        #plotoptimres(mytpd,para,par);
      endfor
      #Output table: prefactor, Ea, coverage, ssq
      row=[v,para(4),para(2),ssq]
      fitret=[fitret; row]
      
      figure(99);
      clf();
      hold on;
      plot(mytpd.T,mytpd.i,'color',mytpd.color);
      plot(mytpd.T,ptotn,'color',mytpd.color);
      plot(mytpd.T,(mytpd.i-ptotn)*10,'color','black');
      plotoptimres(mytpd,para,par);
      #input("Ok for next v?");
      drawnow();
    endfor
    printf("Final values: Ea=%f eV, theta=%f ML\n",para(4),para(2));
    result=fitret;
    save("-text",strcat(mytpd.filename,".fit"),"fitret");
  endif;
endfunction;

function [para,ssq]=optimCoverage(opts,para,mytpd,par)
  fh=@(fp)fitcoverage(fp,para,mytpd,par);
  fp=para(2);
  [paro,ssq]=fminsearch(fh,fp,opts);
  printf("Final SSQ: %e\n",ssq);
  #paro= fminunc(fh,fp,opts);
  para(2)=paro(1);
endfunction;

function [para,ssq]=optimEa(opts,para,mytpd,par)
  printf("Ea optimisation\n");
  fp=para(4);
  fh=@(fp)fitEa(fp,para,mytpd,par);
  [paro,ssq]=fminsearch(fh,fp,opts);
  printf("Final SSQ: %e\n",ssq);
  para(4)=paro(1);
endfunction;

function plotoptimres(mytpd,para,par)
  popt=calcPn(mytpd,para,par);
  plot(mytpd.T,popt,"color","red");
  plot(mytpd.T,(mytpd.i-popt)*10,'color','blue');
  #input("ok, confirm for next point");
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