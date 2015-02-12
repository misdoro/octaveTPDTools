function fitMultiTPD(dats,param)
#Fit with single Ea for different prefactors
#File for energy distribution fit
#Energy distribution parameters - check and display help
#Fit energy distribution for few predefined prefactors
#Fit with various determinded energy distributions, plot sigma(E)

#Check the input data:
filesCount=rows(dats.ordRates);
if (filesCount<2)
  printf("Not enough files to perform multi-rate fit!\n")
  return;
endif;

minrate=dats.ordRates(1,2);
maxrate=dats.ordRates(end,2);
if (minrate*5>maxrate)
  printf("Insufficient range of heating rates!\n");
  return;
endif

#Coarse search of prefactor based on a single Ea:
if (param.debug)
  printf("Coarse 1E search for prefactor\n")
endif

fitpar=loadFitPar();

if (~isfield(fitpar,'estv'))
  figure(1);
  clf;
  vs=getOptionValue(param,"","vs",logspace(10,20,11));
  Efits=[];
  for v=vs
    Eline=[];
    for idx=1:filesCount;
      filename=dats.filenames{dats.ordRates(idx,1)};
      if (param.debug)
        printf("v=%e filename:%s\n",v,filename);      
      endif
    
      mytpd=loadTPD(filename,param);
    
      if (length(mytpd.i)>10)
        fit=initFitParam(mytpd,param,v)
        fiti=fit;
        fit=fitCovScale(mytpd,fit)
        fit=fitE0(mytpd,fit)
        
        if (param.debug)
          clf();
          hold on;
          plot(mytpd.T,mytpd.i,'color','blue');
          plot(mytpd.T,modelTPDmc(mytpd.T,fiti),'color','red');
          plot(mytpd.T,modelTPDmc(mytpd.T,fit),'color','green');
          drawnow();
        endif
        Eline=[Eline,fit.E0]
      else
        printf("No appropriate data in file %s\n",filename);
      endif
    endfor
    Efits=[Efits;Eline]
  endfor
  stds=std(Efits,0,2)
  [minstd,minstdi]=min(stds);
  minv=vs(minstdi)
  fitpar.estv=minv
  save("-text","fit.par","fitpar");
else
  printf("Skipping prefactor coarse search, use saved value of %.1e\n",fitpar.estv);
endif;

#filename=indata.filenames{indata.ordRates(idx,1)};
#printf("\n-----------------\n%s\n",filename);
#load(filename);

endfunction;

function fit=initFitParam(mytpd,param,v)
[maxi,maxii]=max(mytpd.i);
      maxT=mytpd.T(maxii);
      fit.v=v;
      fit.E0=estimEaVTm(v,maxT,mytpd.rate);
      fit.dE=0.1;
      fit.thetas=0.1;
      fit.scale=1;
      fit.ml=2e-8;
      fit.debug=param.debug;
      fit.rate=mytpd.rate;
endfunction

function mytpd=loadTPD(filename,param)
    load(filename);
    smass=param.mass(1);
		mytpd=getMassData(tpd,param.displayT,smass);
    baseline=getOptionValue(param,filename,"bline",0);
    mytpd.i=mytpd.i-baseline;
    factor = round(length(mytpd.i)/100);
    mytpd=decimateTPD(mytpd,factor);
endfunction;

function ret=loadFitPar()
  [f.info, f.err, f.msg]=stat("fit.par");
  if (f.err>=0);
    load "fit.par";
  endif;
  if (exist("fitpar","var"))
    ret=fitpar;
  else
    ret=struct();
  endif
endfunction;
