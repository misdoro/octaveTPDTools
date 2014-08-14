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
  dosecov=doseintg/par.monolayd;

  if (~isfield(par,"mincov")|| ~isfield(par,"maxcov"))
    par.mincov=0.1;
    par.maxcov=0.9;
  endif
  #Check if we have file-specific options defined
  fns=strsplit(mytpd.filename,".");
  fn=fns{1};
  if (isfield(par,fn))
    parfs=getfield(par,fn);
    if (isfield(parfs,"bline"))
      par.bline=parfs.bline;
    endif;
    if (isfield(parfs,"cov"))
      dosecov=parfs.cov;
    endif;
      
    #Cutoff temperatures
    if (isfield(parfs,"minT") && isfield(parfs,"maxT"))
      par.cutT.min=parfs.minT;
      par.cutT.max=parfs.maxT;
    endif
    
    #Cutoff coverage ratio, used by default
    if (isfield(parfs,"mincov")&& isfield(parfs,"maxcov"))
      par.mincov=parfs.mincov;
      par.maxcov=parfs.maxcov;
    endif;
      
  endif
    
  if (isfield(par,"bline"))
    mytpd.i=mytpd.i-par.bline;
  endif;
  if (isfield(par,"cutT"))
    mytpd=cutTemp(mytpd,par.cutT);
  else
    cov=cumtrapz(mytpd.t,mytpd.i);
    covnorm=cov./cov(end);
    cutstart=max(find(covnorm<par.mincov));
    cutend=min(find(covnorm>par.maxcov));
    mytpd=cutIndex(mytpd,length(mytpd.i),cutstart,cutend);
  endif;
  #Decimate TPD to reduce to defined number of points, if asked
  if (isfield(par,"decimate") && par.decimate)
    factor = round(length(mytpd.i)/par.np)
    mytpd=decimateTPD(mytpd,factor);
  endif
	
	para=[par.v,dosecov,par.monolay*mytpd.rate,par.E];
	ptotn=calcPn(mytpd,para,par);
	
  
	
  figure(param.fig.model);
  plot(mytpd.T,ptotn,'color',mytpd.color,'linewidth',2);
	plot(mytpd.T,(mytpd.i-ptotn)*10,'color',mytpd.color);
  plot(mytpd.T,mytpd.i,".",'color',mytpd.color);
	drawnow();
  
  ssq=fitcoverage((dosecov),para,mytpd,par);
  printf("Initial SSQ: %e\n",ssq);
  
  
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
        
        
        #opts=optimset("Display","iter","TolX",1e-9);
        #opts=optimset("Display","final","TolX",1e-9);
        opts=optimset("Display","final","TolX",1e-5);
        
        [para,ssq]=optimEa(opts,para,mytpd,par);
        #plotoptimres(mytpd,para,par);
        
        [para,ssq]=optimCoverage(opts,para,mytpd,par);
        #plotoptimres(mytpd,para,par);
      endfor
      #Output table: prefactor, Ea, coverage, ssq
      row=[v,para(4),para(2),ssq]
      fitret=[fitret; row]
      
      figure(param.fig.modelfit);
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
    #Output table: prefactor, Ea, coverage, ssq
    save("-text",strcat(mytpd.filename,".fit"),"fitret");
  endif;
endfunction;

function [para,ssq]=optimCoverage(opts,para,mytpd,par)
  printf("Coverage optimisation\n");
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
  if (fp(1)<0)
    fp(1)=eps;
  endif;
  parv(2)=fp(1)
	p=calcPn(tpd,parv,par);
	ssq=sumsq(p-tpd.i);
endfunction;

function ssq=fitEa(fp,parv,tpd,par)
  if (fp(1)<0)
    fp(1)=eps;
  endif;
  parv(4)=fp(1)
	p=calcPn(tpd,parv,par);
	ssq=sumsq(p-tpd.i);
endfunction;