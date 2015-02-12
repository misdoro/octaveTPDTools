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
        Eline=[Eline,fit.E0];
      else
        printf("No appropriate data in file %s\n",filename);
      endif
    endfor
    Efits=[Efits;Eline]
  endfor
  stds=std(Efits,0,2);
  [minstd,minstdi]=min(stds);
  minv=vs(minstdi);
  fitpar.estv=minv;
  printf("Found optimum prefactor to be %.1e\n",fitpar.estv);
  save("-text","fit.par","fitpar");
else
  printf("Skipping prefactor coarse search, use saved value of %.1e\n",fitpar.estv);
endif;

#Fit ea distribution for Vopt, Vopt*100, Vopt/100
if (param.debug)
  printf("Fit Ea distributions for chosen prefactor values\n")
endif
if (~isfield(fitpar,'defits'))
  defits={}
  figure(1)
  idx=0;
  for v=[fitpar.estv/100,fitpar.estv,fitpar.estv*100]
    #Work with the lowest-rate TPD to get the most accurate data.
    filename=dats.filenames{dats.ordRates(1,1)};
    mytpd=loadTPD(filename,param);
    fit=initFitParam(mytpd,param,v);
    fit.E0-=4*fit.dE;
    fit.thetas=0.01*ones(15,1);
    fiti=fit
    fit=fitPartCoverages(mytpd,fit)
    
    defits{++idx}=fit;
    clf;
    hold on;
    p=modelTPDmc(mytpd.T,fiti);
    po=modelTPDmc(mytpd.T,fit);
    plot(mytpd.T,p,'color','red');
    plot(mytpd.T,mytpd.i,'color','blue');
    plot(mytpd.T,po,'color','green');
    drawnow();
  endfor
  fitpar.defits=defits;
  save("-text","fit.par","fitpar");
else
  printf("Using predefined energy distributions\n");
  fitpar.defits
endif

#Fine fit of prefactor using energy distribution
figure(2);
clf;
hold on;
colors={'green','blue','red'}
if (~isfield(fitpar,'stds'))
vs=logspace(log10(fitpar.estv/1000),log10(fitpar.estv*1000),11)

for vidx=[1,2,3]
  Efits=[];
  for v=vs
    Eline=[];
    for idx=1:filesCount;
      filename=dats.filenames{dats.ordRates(idx,1)};
      if (param.debug)
        printf("step vidx=%d ",vidx);
        printf("v=%.1e filename:%s\n",v,filename);      
      endif 
      mytpd=loadTPD(filename,param);
      fit=fitpar.defits{vidx}
      fit.v=v;
      fit.rate=mytpd.rate;
      fiti=fit;
      fit=fitE0(mytpd,fit);
      fit=fitCovScale(mytpd,fit);
      fit=fitE0(mytpd,fit);
      fit=fitCovScale(mytpd,fit);
      if (param.debug)
      figure(1);
        clf();
        hold on;
        plot(mytpd.T,mytpd.i,'color','blue');
        plot(mytpd.T,modelTPDmc(mytpd.T,fiti),'color','red');
        plot(mytpd.T,modelTPDmc(mytpd.T,fit),'color','green');
        drawnow();
      endif
      Eline=[Eline,fit.E0];
    endfor;
    Efits=[Efits;Eline]
    
  endfor;
  stds{vidx}=std(Efits,0,2);
  legendtxt{vidx}=sprintf('v=%.1e',fitpar.defits{vidx}.v);
  figure(2);
  hold on;
  semilogx(vs,stds{vidx},'color',colors{vidx});
  legend(legendtxt);
  drawnow;
  fitpar.stds=stds;
  save("-text","fit.par","fitpar");
endfor

else
  printf("Using previous fit errors\n");
  fitpar.stds
  colors={'green','blue','red'}
  figure(2);
  clf();
  hold on;
  vs=logspace(log10(fitpar.estv/1000),log10(fitpar.estv*1000),11)
  
  for vidx=[1,2,3]
    legendtxt{vidx}=sprintf('v=%.1e',fitpar.defits{vidx}.v);
    semilogx(vs,fitpar.stds{vidx},'color',colors{vidx});
    legend(legendtxt);
    drawnow;
  endfor;
endif

endfunction;

function fit=initFitParam(mytpd,param,v)
  [maxi,maxii]=max(mytpd.i);
  maxT=mytpd.T(maxii);
  fit.v=v;
  fit.E0=estimEaVTm(v,maxT,mytpd.rate);
  fit.dE=0.01;
  fit.thetas=0.1;
  fit.scale=1;
  fit.ml=2e-8;
  fit.debug=param.debug;
  fit.rate=mytpd.rate;
endfunction

function mytpd=loadTPD(filename,param,decimate=1)
  load(filename);
  smass=param.mass(1);
  mytpd=getMassData(tpd,param.displayT,smass);
  baseline=getOptionValue(param,filename,"bline",0);
  mytpd.i=mytpd.i-baseline;
  mytpd.i=supsmu(mytpd.T,mytpd.i,'spa',0.05);
  if (decimate)
    factor = round(length(mytpd.i)/50);
    mytpd=decimateTPD(mytpd,factor);
  endif;
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
