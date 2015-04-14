function fitMultiTPD(dats,param)
#Desorption prefactor search algorithm:
#-Fit with single Ea for different prefactors
#-Fit energy distribution for few predefined prefactors, based on the TPD with lowest rate
#-Fit varying prefactor with 3 determinded energy distributions, plot sigma(E)

#Check the input data:
filesCount=rows(dats.ordRates);
if (filesCount<2)
  printf("Not enough files to perform multi-rate fit!\n")
  return;
endif;

minrate=dats.ordRates(1,2);
maxrate=dats.ordRates(end,2);
if (minrate*3>=maxrate)
  printf("Insufficient range of heating rates!\n");
  return;
endif

#Coarse search of prefactor based on a single Ea:
if (param.debug)
  printf("Coarse 1E search for prefactor\n")
endif

fitpar=loadFitFile();

if (~isfield(fitpar,'estv'))
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
          figure(getFigIndex("fit_estv"));
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
  
  figure(getFigIndex("fit_estv_stds"));
  hold on;
  semilogx(vs,stds,'color','black');
  xlabel("Prefactor value");
  ylabel("Fitted energy stddev");
  
  [minstd,minstdi]=min(stds);
  minv=vs(minstdi);
  fitpar.estv=minv;
  printf("Found optimum prefactor to be %.1e\n",fitpar.estv);
  save("-text","fit.par","fitpar");
else
  printf("Skipping prefactor coarse search, use saved value of %.1e\n",fitpar.estv);
  figure(getFigIndex("fit_1Efit"));
  hold on;
  xlabel("Temperature, K");
  ylabel("desorption signal, ML/K");
  for idx=1:filesCount;
      filename=dats.filenames{dats.ordRates(idx,1)};
      mytpd=loadTPD(filename,param);
      fit=initFitParam(mytpd,param,fitpar.estv)
      fit=fitCovScale(mytpd,fit)
      fit=fitE0(mytpd,fit)
      legendtxt{idx*2-1}=sprintf("exp, %d K/min",round(mytpd.rate*60));
      legendtxt(idx*2)=sprintf("mod. Eo=%.3f,v=%.1e, %d K/min",fit.E0,fitpar.estv,round(mytpd.rate*60));
      ydat=mytpd.i./(param.monolayer.*(mytpd.rate));
      plot(mytpd.T,ydat,'color','blue',"marker","x");
      ydat=modelTPDmc(mytpd.T,fit)./(param.monolayer.*(mytpd.rate));
      plot(mytpd.T,ydat,'color','green');
      drawnow();
  endfor;
  legendtxt
  h=legend(legendtxt);
  set (h, 'fontsize', 10);
endif;

#Fit ea distribution for Vopt, Vopt*100, Vopt/100
if (param.debug)
  printf("Fit Ea distributions for chosen prefactor values\n")
endif
if (~isfield(fitpar,'defits'))
  defits={}
  idx=0;
  for v=[fitpar.estv/100,fitpar.estv,fitpar.estv*100]
    [fit,fiti,mytpd]=fitEnergyDistribution(dats,param,v);
    
    defits{++idx}=fit;
    if (param.debug)
      figure(getFigIndex("fit_debug"));
      clf;
      hold on;
      p=modelTPDmc(mytpd.T,fiti);
      po=modelTPDmc(mytpd.T,fit);
      plot(mytpd.T,p,'color','red');
      plot(mytpd.T,mytpd.i,'color','blue');
      plot(mytpd.T,po,'color','green');
      drawnow();
      #input("is the fit OK?");
    endif
  endfor
  fitpar.defits=defits;
  save("-text","fit.par","fitpar");
else
  printf("Using predefined energy distributions\n");
  fitpar.defits
endif

#Fine fit of prefactor using energy distribution
colors={'green','blue','red'}
if (~isfield(fitpar,'stds'))
  vs=logspace(log10(fitpar.estv/100),log10(fitpar.estv*100),20)
  fitpar.vsfine=vs;
  legendtxt{1}="1E fit";
  for vidx=[1,2,3]#Iterate over energy distribution fits
   Efits=[];
   
   for v=vs#Iterate over prefactors in range
     Eline=[];
     for idx=1:filesCount;#Iterate over available files
       filename=dats.filenames{dats.ordRates(idx,1)};
       if (param.debug)
         printf("step vidx=%d ",vidx);
         printf("v=%.1e filename:%s\n",v,filename);      
       endif 
       mytpd=loadTPD(filename,param);
       fit=fitpar.defits{vidx}
       fit.v=v;
       fit.E0=estimE0(mytpd,fit)-4*fit.dE;
       fit.rate=mytpd.rate;
       fiti=fit;
       fit=fitE0(mytpd,fit);
       fit=fitCovScale(mytpd,fit);
       fit=fitE0(mytpd,fit);
       fit=fitCovScale(mytpd,fit);
       if (param.debug)
         figure(getFigIndex("fit_debug"));
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
   fitpar.stds{vidx}=std(Efits,0,2);
   
   legendtxt{vidx+1}=sprintf('v=%.1e',fitpar.defits{vidx}.v);
   figure(getFigIndex("fit_estv_stds"));
   clf();
   hold on;
   xlabel("Prefactor value");
   ylabel("Fitted energy stddev");
   semilogx(vs,fitpar.stds{vidx},'color',colors{vidx});
   h=legend(legendtxt);
   set (h, 'fontsize', 10);
   drawnow;
   endfor
else
  printf("Using previous fit errors\n");
  fitpar.stds
endif
  colors={'green','blue','red'}
  figure(getFigIndex("fit_estv_stds"));
  clf();
  hold on;
  vs=fitpar.vsfine;
  mins=[];
  legendtxt={}
  for vidx=[1,2,3]
    legendtxt{vidx}=sprintf('v=%.1e',fitpar.defits{vidx}.v);
    semilogx(vs,fitpar.stds{vidx},'color',colors{vidx});
    
    [minstd,minstdi]=min(fitpar.stds{vidx});
    minv=vs(minstdi);
    mins=[mins,minv];
  endfor;
  h=legend(legendtxt);
  set (h, 'fontsize', 10);
  ylabel("Energy fit stddev");
  xlabel("Prefactor value");
  drawnow;
  fitpar.bestv=mean(mins)
  save("-text","fit.par","fitpar");


if (~isfield(fitpar,'fitdE')||isfield(param,'refitdE'))
  #Final fit of the energy distribution
  printf("Fitting energy distribution for the best prefactor of v=%.1e\n",fitpar.bestv);
  [fitpar.fitdE,fiti,mytpd]=fitEnergyDistribution(dats,param,fitpar.bestv);
  save("-text","fit.par","fitpar");
  figure(getFigIndex("fit_bestfit"));
  clf;
  hold on;
  p=modelTPDmc(mytpd.T,fiti);
  po=modelTPDmc(mytpd.T,fitpar.fitdE);
  plot(mytpd.T,p,'color','red');
  plot(mytpd.T,mytpd.i,'color','blue');
  plot(mytpd.T,po,'color','green');
  drawnow();
endif
  printf("Showing previously prepared dE fit\n");
  fit=fitpar.fitdE
  figure(getFigIndex("fit_finaldE"));
  clf();
  hold on;
  np=length(fit.thetas);
  Epts=linspace(fit.E0,fit.E0+(np-1)*fit.dE,np);
  bar(Epts,fit.thetas,"linewidth",2);
  xlabel("Ea, eV");
  ylabel("Theta_i, ML");

  figure(getFigIndex("fit_dEfit"));
  hold on;
  xlabel("Temperature, K");
  ylabel("desorption signal, ML/K");
  for idx=1:filesCount;
      filename=dats.filenames{dats.ordRates(idx,1)};
      mytpd=loadTPD(filename,param);
      fit=fitpar.fitdE;
      fit.rate=mytpd.rate;
      fit=fitCovScale(mytpd,fit);
      legendtxt{idx*2-1}=sprintf("exp, %d K/min",round(mytpd.rate*60));
      legendtxt(idx*2)=sprintf("mod.1E,v=%.1e, %d K/min",fit.v,round(fit.rate*60));
      ydat=mytpd.i./(param.monolayer.*(mytpd.rate));
      plot(mytpd.T,ydat,'color','blue',"marker","x");
      ydat=modelTPDmc(mytpd.T,fit)./(param.monolayer.*(mytpd.rate));
      plot(mytpd.T,ydat,'color','green');
      drawnow();
  endfor;
  h=legend(legendtxt);
  set (h, 'fontsize', 10);
  
endfunction;

function [fit,fiti,mytpd]=fitEnergyDistribution(dats,param,v)
  #Work with the lowest-rate TPD to get the most accurate data.
  filename=dats.filenames{dats.ordRates(1,1)};
  mytpd=loadTPD(filename,param);
  fit=initFitParam(mytpd,param,v);
  shift=6;
  fit.E0-=shift*fit.dE;
  fit.thetas=0.01*ones(fit.np+shift,1);
  fiti=fit
  fit=fitPartCoverages(mytpd,fit)
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
