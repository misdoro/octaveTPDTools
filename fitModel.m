function result=fitModel(mytpd,param,result,press,dose,fitparam);
	#par.mids=mytpd.mids;
	
	#doseintg=0;
	#if (isfield(mytpd,"version") && mytpd.version>=20140120 && ~mytpd.model)
	#	doseintg=calculateDoseIntegral(getMassData(dose,[],param.selectedmass));
	#else
	#	doseintg=mytpd.doseintg;
	#endif

	#par=loadParamFile(par);
  #dosecov=doseintg/par.monolayd;
    
  #if (isfield(par,"bline"))
  #  mytpd.i=mytpd.i-par.bline;
  #endif;
  
  if (isfield(par,"cutT"))
    mytpd=cutTemp(mytpd,par.cutT);
  else
    cov=cumtrapz(mytpd.t,mytpd.i);
    covnorm=cov./cov(end);
    mincov=getOptionValue(par,mytpd.filename,"mincov",0.1);
    maxcov=getOptionValue(par,mytpd.filename,"maxcov",0.9);
    cutstart=max(find(covnorm<mincov));
    cutend=min(find(covnorm>maxcov));
    mytpd=cutIndex(mytpd,length(mytpd.i),cutstart,cutend);
  endif;
  
  #Decimate TPD to reduce to defined number of points, if asked
  #if (isfield(par,"decimate") && par.decimate)
  #  factor = round(length(mytpd.i)/par.np)
  #  mytpd=decimateTPD(mytpd,factor);
  #endif
	
	#para=[par.v,dosecov,par.monolay*mytpd.rate,par.E,par.Es];
	#ptotn=calcPn(mytpd,para,par);
	
  
	
  #figure(param.fig.model);
  #plot(mytpd.T,ptotn,'color',mytpd.color,'linewidth',2);
	#plot(mytpd.T,(mytpd.i-ptotn)*10,'color',mytpd.color);
  #plot(mytpd.T,mytpd.i,".",'color',mytpd.color);
	#drawnow();
  
  #ssq=fitcoverage((dosecov),para,mytpd,par);
  #printf("Initial SSQ: %e\n",ssq);
  
  
	if (index(param.tools,'F'))
    if isfield(par,"vs")
      vs=par.vs
    else
      vs=par.v
    endif
    fitret=[];
    for v=vs
      fit.v=v;
      #para(1)=v
      #printf("Iteration for v=%e\n",v);
      #ssq=0;
      for i=[1,2]
        
        
        #opts=optimset("Display","iter","TolX",1e-9);
        #opts=optimset("Display","final","TolX",1e-9);
        opts=optimset("Display","final","TolX",1e-5);
        
        #[para,ssq]=optimEa(opts,para,mytpd,par);
        #plotoptimres(mytpd,para,par);
        
        #[para,ssq]=optimCoverage(opts,para,mytpd,par);
        #plotoptimres(mytpd,para,par);
        
        if (isfield(par,"fitde")&& par.fitde~=0)
          [para,ssq]=optimEs(opts,para,mytpd,par);
        endif
      endfor
      #Output table: prefactor, Ea, coverage, ssq
      row=[v,para(4),para(2),ssq,para(5)]
      fitret=[fitret; row]
      
      #figure(param.fig.modelfit);
      #clf();
      #hold on;
      #plot(mytpd.T,mytpd.i,'color',mytpd.color);
      #plot(mytpd.T,ptotn,'color',mytpd.color);
      #plot(mytpd.T,(mytpd.i-ptotn)*10,'color','black');
      #plotoptimres(mytpd,para,par);
      #input("Ok for next v?");
      #drawnow();
    endfor
    #printf("Final values: Ea=%f eV, theta=%f ML\n",para(4),para(2));
    result=fitret;
    #Output table: prefactor, Ea, coverage, ssq, Es
    save("-text",strcat(mytpd.filename,".fit"),"fitret");
  endif;
endfunction;