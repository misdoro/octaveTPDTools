#!/usr/bin/octave -qi
#Arguments: mid, min T, max T, display, monolayer

if (nargin==0)
	printf("\n\
  ################################################################\n\
  # TPD data treatment program, LPMAA/LERMA, UPMC, Paris         #\n\
  ################################################################\n\
  #Usage: display.m act [M1,M2,MN startTemp endTemp monolayer]   #\n\
  #Actions:                                                      #\n\
  #-----------display --------------                             #\n\
  #d = plot TPDs                                                 #\n\
  #r = plot real non-smooth data                                 #\n\
  #x = extract the baseline from the data                        #\n\
  #c = use coverage-based colors                                 #\n\
  #N = Normalize to ML/K as y axis for display                   #\n\
  #D = plot dose graph                                           #\n\
  #f = plot the temperature fit quality                          #\n\
  #p = plot pressure                                             #\n\
  #P = plot pressure vs Iqms                                     #\n\
  #R = plot the IR spectrum, add points to the TPD plot          #\n\
  #-----------Ea extract -----------                             #\n\
  #l = log plot of i over 1/T                                    #\n\
  #e = energy estimation using inversion plot over coverage      #\n\
  #t = plot temperature points on the inversion curves every 5K  #\n\
  #V = perform a search of prefactor based on multi-rate TPDs    #\n\
  #E = Fit and plot initial coverages for all available TPDs     #\n\
  #m = model TPD with 1-st order process                         #\n\
  #F = Fit model parameters                                      #\n\
  #-----------Fine tune, misc.------                             #\n\
  #I = interactive parameters                                    #\n\
  #i = print files info: available masses and T range            #\n\
  #v = ask to clear plot after each iteration                    #\n\
  #T = treat isotherm desorption                                 #\n\
  #u = run the user init.m, iterate.m and final.m, if available  #\n\
  #s = save the graph points in ascii format                     #\n\
  ################################################################\n\n");
endif

param.monolayer=2e-08;#Xe /amorph 1,2e-6 A*s 100K
param.monolayer=useArgument(argv(),5,param.monolayer);
param.debug=0;
param=loadParamFile(param);

datindex.filenames=findDatFiles(pwd);
#By default, use the last MID (most precise timestamp) and available 
#temperature range from the first input file.
param=getFileInfo(datindex,param);
param.displayT.min=useArgument(argv(),3,param.displayT.min);
param.displayT.max=useArgument(argv(),4,param.displayT.max);
param.mass=useArgument(argv(),2,param.mass);
for midx=1:length(param.mass)
	printf("Displayed temperature: %.2f to %.2f K\nDisplayed MID: %d\n"
	,param.displayT.min,param.displayT.max,param.mass(midx));
end

datindex=indexDatFiles(datindex,param);

param.tools=useArgument(argv(),1,"i");

pkg load optim;

#################################################
# Plot TPD
#################################################
function result=plotTPD(mytpd,param,result,press,dose);
	mini=min(mytpd.i);
	if (index(param.tools,'r'))
		mytpd.i_sm=mytpd.i;
	else
    smooth=getOptionValue(param,mytpd.filename,"smooth",0.005);
		mytpd.i_sm=supsmu(mytpd.T,mytpd.i,'spa',smooth);
	endif;
  #default plot style
	ls="-";
  lw=2;
  ms="none";
  
	if (isfield(mytpd,"model")&& mytpd.model>0);
    #Dotted line for model TPD
		ls=":";
	endif;
  if (isfield(param,"points"))
    #Decimated points from par.points=N
	  ls="none";
    ms=".";
    lw=1;
    ps=param.points
	  plt=plot(mytpd.T(1:ps:end),mytpd.i_sm(1:ps:end),"linewidth",lw,"linestyle",ls,"marker",ms,"color",mytpd.color);
  else  
    #Normal plot
    plt=plot(mytpd.T,mytpd.i_sm,"linewidth",lw,"linestyle",ls,"marker",ms,"color",mytpd.color);
  endif
	
	[maxi,maxidx]=max(mytpd.i_sm);
	maxT=mytpd.T(maxidx);
	if (length(param.mass)>1)
		txt=sprintf("%d:M=",mytpd.idx,mytpd.mass)
  else
    txt=sprintf("%d",mytpd.fileidx)
	endif
	text(maxT,maxi,txt);
	fn=strrep(mytpd.filename,"_","-");
	if (isfield(param,"publish"))
    legendtext=sprintf("%s: %3.2f ML @ %d K/min ",txt,mytpd.intg/param.monolayer,round(mytpd.rate*60));
  else
    legendtext=sprintf("%s:%s(%3.2f ML)",txt,fn,mytpd.intg/param.monolayer);
  endif
	result=retAppend(result,"legend",legendtext);
  result=retAppend(result,"plothdl",plt);
	
	doseintg=0;
	if (isfield(mytpd,"version") && mytpd.version>=20140120 && ~mytpd.model)
		dosdat=getMassData(dose,[],param.selectedmass);
		doseintg=calculateDoseIntegral(dosdat,0);
	else
		doseintg=mytpd.doseintg;
	endif
	result=retAppend(result,"doses",doseintg);
	result=retAppend(result,"integrals",mytpd.intg);
  if (mytpd.baseline>0)
    printf("Corrected ")
  endif
	printf("TPD integral: %.3e\n",mytpd.intg);
	printf("Dose integral: %.3e\n",doseintg);
	printf("Itpd/Idose: %.3f\n",mytpd.intg/doseintg);
  
  #Rate and maximum position
  printf("TPD rate %.2f K/s (%.1f K/min)\n",mytpd.rate,mytpd.rate*60);
  [maxi, maxin]=max(mytpd.i_sm);
  printf("Temperature of maximum: %3.2f K\n",mytpd.T(maxin));
endfunction

if (index(param.tools,'d'))
	figure(getFigIndex("disp"));
	hold on;
	ret=iterateTpd(datindex,param,@plotTPD);
  if (index(param.tools,'N'))
    ylabel("Desorption flow (ML/K)");
  else
	  ylabel("Desorption flow (arb.u.)");
  endif
	xlabel("Temperature (K)");
  
	if (isfield(ret,"legend"))
    legloc="northeast";
    if (isfield(param,"legendloc"))
      legloc=param.legendloc;
    endif;
		#legend("boxon","right");
		h=legend(ret.legend,"location",legloc);
		set (h, 'fontsize', 10);
	endif;
	xlim([param.displayT.min,param.displayT.max]);
  lim=ylim();
  lim(1)=-0.025*abs(lim(2));
  ylim(lim);
	
	#Plot doses
	if (isfield(ret,"doses"))
		figure(getFigIndex("doses"));
		plotDoses(ret);
	endif;
endif;

#################################################
# iterate with user-provided code
#################################################
function result=iterTPD(mytpd,param,result,press,dose);
  source("iterate.m")
  result=result;
endfunction;

if (index(param.tools,'u'))
  [init.info, init.err, init.msg]=stat("init.m");
  [iter.info, iter.err, iter.msg]=stat("iterate.m");
  [fin.info, fin.err, fin.msg]=stat("final.m");
	if (iter.err>=0);
    #Prepare the user figure
    figure(getFigIndex("user"));
	
    #First check if we have the iteration script
    if (init.err>=0);
      source("init.m");
    endif;
    ret=iterateTpd(datindex,param,@iterTPD);
    
    if (fin.err>=0);
      source("final.m");
    endif;
  endif;
endif;

#################################################
# Plot extended dose information
#################################################
function result=plotDoses(mytpd,param,result,press,dose);
	dose=getMassData(dose,[],param.selectedmass);
  normdosi=dose.i*1e10;
  figure(getFigIndex("doseext"));
	semilogy(dose.t,normdosi,"linewidth",2,"color",mytpd.color);
	[dosi,endpt]=calculateDoseIntegral(dose,0);
  plot(dose.t(endpt),normdosi(endpt),"o")
	printf("Dose temperature %.2f K\n", mean(dose.T));
  printf("Dose integral %.2e\n",dosi);
	
	[maxi,maxidx]=max(normdosi);
	maxt=dose.t(maxidx);
	#text(maxt,maxi,num2str(mytpd.fileidx));
  
endfunction;


if (index(param.tools,'D'))
	figure(getFigIndex("doseext"));
	hold on;
	ylabel("Dose current, Ax10^-10");
	xlabel("Time (s)");
	ret=iterateTpd(datindex,param,@plotDoses);
endif;
#################################################
# Plot gauge pressure
#################################################
function result=plotP(mytpd,param,result);
	
	plot(mytpd.T,mytpd.pi,"linewidth",2,"color",mytpd.color);
	
	[maxi,maxidx]=max(mytpd.pi);
	maxT=mytpd.T(maxidx);
	text(maxT,maxi,strcat(num2str(mytpd.idx)));
endfunction

if (index(param.tools,'p'))
	figure(getFigIndex("press"));
	hold on;
	ret=iterateTpd(datindex,param,@plotP);
	ylabel("Pressure (torr)");
	xlabel("Temperature (K)");
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	
endif;

#################################################
# Plot QMS current vs pressure 
#################################################
if (index(param.tools,'P'))
	figure(getFigIndex("qipress"));
	hold off;
	ylabel("Iqms, A");
	xlabel("Pressure (torr)");
	ret=iterateTpd(datindex,param,@dispplotqP);
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	
endif;

#################################################
# Plot temperature fit quality
#################################################
if (index(param.tools,'f'))
	figure(getFigIndex("Tfitq"));
	hold on;
	ylabel("T ramp fit error, K");
	xlabel("T, K");
  
  figure(getFigIndex("Tramp"));
	hold on;
	ylabel("Temperature, K");
	xlabel("Time, s");
  
	ret=iterateTpd(datindex,param,@dispTfitq);
  
	figure(getFigIndex("Tfitq"));
  if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
  
  figure(getFigIndex("Tramp"));
  if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
endif;

##################################################
# Log(i) over 1/T plot + linear Eads extraction  #
##################################################
function result=plotInvT(mytpd,param,result);
	mytpd.invT=1000./mytpd.T;
	mytpd.logi=log(mytpd.i);
	
	figure(getFigIndex("log"));
	plot(mytpd.invT,mytpd.logi,"linewidth",2,"color",mytpd.color);

	[eads,lnv,win]=findLogEAds(mytpd);
  eads=eads*1000;%Correct because we use 1000/T
	eVKcalmul=23.0609;
  eVKJmul=96.4869;
	printf("Eads = %f eV, %f kcal/mol, %f kJ/mol\n",eads,eads*eVKcalmul,eads*eVKJmul);
	txtx=mytpd.invT(1);
	txty=mytpd.logi(1);
	text(txtx,txty,strcat("<",num2str(mytpd.idx)));
	
	if(getFigIndex("disp",0))
		figure(getFigIndex("disp"));
		plot(mytpd.T(win.istart),mytpd.i(win.istart),"o");
		plot(mytpd.T(win.iend),mytpd.i(win.iend),"o");
	endif;
	
  ltxt=sprintf("exp. θ_0=%.2f ML",mytpd.intg/param.monolayer);
  result=retAppend(result,"legend",ltxt);
  ltxt=sprintf("lin. fit,Ea=%.2f eV",eads);
  result=retAppend(result,"legend",ltxt);
  ltxt="fit error x1000";
  result=retAppend(result,"legend",ltxt);
  
  
endfunction

if (index(param.tools,'l'))
	figure(getFigIndex("log"));
	hold on;
	ylabel("ln(Φdes)")
	xlabel("Inverse Temperature (1000/T)")
	
	ret=iterateTpd(datindex,param,@plotInvT);
	
	if (isfield(ret,"legend"))
    figure(getFigIndex("log"));
		legend(ret.legend);
	endif;
	
endif;


##################################################
#Eads inversion	                          			 #
##################################################
if (index(param.tools,'e'))
	figure(getFigIndex("eads"));
	hold on;
	ylabel("Ea (eV)")
	xlabel("Coverage (ML)")
	ret3=iterateTpd(datindex,param,@plotEAds);
endif

#################################################################
# Fit multi-rate TPDs to find prefactor and energy distribution #
#################################################################
if (index(param.tools,'V'))
  fitMultiTPD(datindex,param);
endif



#######################################################
# Fit initial coverages for known prefactor           #
#######################################################
if (index(param.tools,'E'))
  figure(getFigIndex("covsites"));
	hold on;
  figure(getFigIndex("covfiterr"));
	hold on;
  result=iterateTpd(datindex,param,@fitMultiCoverages);
  if (isfield(result,"fitemaxde"))
    printf("Adsorption energies summary\n");
    printf("  θ, ML  ;  Ea, eV;   dEa, eV\n");
    printf("%3.2f ; %5.3f ; %5.3f \n",result.fitemaxde');
  endif
  figure(getFigIndex("covsites"));
  xlabel("Ea, eV");
  ylabel("Sites population");
  figure(getFigIndex("covfiterr"));
  xlabel("T,K");
  ylabel("Desorption fit error, %max");
endif

#######################################################
# TPD modeling, uses (user-defined) function fitModel #
#######################################################

if (index(param.tools,'m'));
	figure(getFigIndex("model"));
	hold on;
	pkg load odepkg;
	iterateTpd(datindex,param,@fitModel);
  figure(getFigIndex("model"));
	xlabel("Temperature (K)")
	ylabel("Desorption signal");
  
  #Plot fit results, if any
  fitfiles=findDatFiles(pwd,"fit");
  if (length(fitfiles)>1)
    figure(getFigIndex("modelediff"));
    toplotv=[];
    toplote=[];
    for filen=fitfiles;
      filen{1}
      load(filen{1});
      if (length(fitret(:,1)==length(toplotv))||length(toplotv)==0)
        toplotv=fitret(:,1);
        toplote=[toplote,fitret(:,2)];
      endif
    endfor
    toplotstd=std(toplote,0,2);
    semilogx(toplotv,toplotstd,"linewidth",2);
    xlabel("Prefactor v");
    ylabel("E stddev(ramp)");
  endif
endif

########################################################
# Print file info for all dat files in the current dir #
########################################################

if (index(param.tools,'i'));
	printInfo(datindex,param);
endif

########################################################
# Experimental isotherm treatment, not complete        #
########################################################
	
if (index(param.tools,'T'));
	figure(getFigIndex("isotherm"));
	hold on;
	baseparam=iterateTpd(datindex,param,@treatIsotherm);
endif

########################################################
# Infra-red spectrum plot, with colour points on TPD   #
########################################################
if (index(param.tools,'R'));
  figure(getFigIndex("IR"));
  ylabel("Absorbance")
	xlabel("Wavenumber (cm-1)")
  set (gca (), "xdir", "reverse");
	hold on;
  [f.info, f.err, f.msg]=stat("IR.irdat");
	if (f.err>=0);
    load ("IR.irdat");
	  points=iterateTpd(datindex,param,@plotIRData,IRDAT);
  endif
endif

########################################################
# Apply pressure correction                            #
########################################################
function result=extractBaseLine(mytpd,param,result,press,dose);
	result=plotTPD(mytpd,param,result,press,dose);
	mytpd.baseLine=interpBaseLine(param.baseLine,mytpd,press);
	mytpd.i=mytpd.i-mytpd.baseLine;
	printf("Corrected TPD integral %.3e\n",trapz(mytpd.t, mytpd.i));
	mytpd.color='black';
	result=plotTPD(mytpd,param,result,press,dose);
	result.Ts{mytpd.idx}=mytpd.T;
	result.is{mytpd.idx}=mytpd.i;
endfunction;


#Save all figures in the end, since they may be updated during the processing
if(getFigIndex("nonsense")>1)
  input("Adjust figures if needed, then press enter to save, ctrl-c to exit");
  saveAsc=index(param.tools,'s')
  format="png";
  if (isfield(param,"imgformat"))
    format=param.imgformat;
  endif;
    
  saveFig("disp","desorption",saveAsc,format);
  saveFig('doses','sticking',saveAsc,format);
  saveFig('doseext','doses',saveAsc,format);
  saveFig('log','logplot',saveAsc,format);
  saveFig('eads','eadsest',saveAsc,format);
  saveFig('press','pressure',saveAsc,format);
  saveFig('qipress','qms-pressure',saveAsc,format);
  saveFig('user','user',saveAsc,format);
  saveFig('IR','FTIR',saveAsc,format);
  saveFig('modelediff','prefactor',saveAsc,format);
  saveFig('Tfitq','Trampqual',saveAsc,format);
  saveFig('Tramp','Tramp',saveAsc,format);
  
  saveFig('fit_1Efit',"VSearch1Efit",saveAsc,format);
  saveFig('fit_estv_stds',"VSearch1Estd",saveAsc,format);
  saveFig('fit_dEfit',"VSearchdEfit",saveAsc,format);
  saveFig('fit_estv_stdsde',"VSearchEstd",saveAsc,format);
  saveFig('fit_finaldE',"VSearchdE",saveAsc,format);
  
  saveFig('covsites',"dEFitSites",saveAsc,format);
  saveFig('covfiterr',"dEFitErr",saveAsc,format);
endif;
exit(0);

