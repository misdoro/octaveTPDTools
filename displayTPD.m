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
param.figindex=0;%Current figure index, to be autoincremented by processings

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
	ls="-";
	if (isfield(mytpd,"model")&& mytpd.model>0);
		ls=":";
	endif;
	
	plot(mytpd.T,mytpd.i_sm,"linewidth",2,"linestyle",ls,"color",mytpd.color);
	
	[maxi,maxidx]=max(mytpd.i_sm);
	maxT=mytpd.T(maxidx);
	txt=num2str(mytpd.idx);
	if (length(param.mass)>1)
		txt=strcat(txt,":M=",num2str(mytpd.mass))
	endif;
	text(maxT,maxi,txt);
	fn=strrep(mytpd.filename,"_","-");
	#legendtext=strcat(txt,":",fn,"(",num2str(mytpd.intg/param.monolayer,"%3.2f"),"ML)");
	legendtext=sprintf("%d: %d K/min (%3.2f ML)",mytpd.idx,mytpd.rate*60,mytpd.intg/param.monolayer);
	result=retAppend(result,"legend",legendtext);
	
	doseintg=0;
	if (isfield(mytpd,"version") && mytpd.version>=20140120 && ~mytpd.model)
		dosdat=getMassData(dose,[],param.selectedmass);
		doseintg=calculateDoseIntegral(dosdat,0);
	else
		doseintg=mytpd.doseintg;
	endif
	result=retAppend(result,"doses",doseintg);
	result=retAppend(result,"integrals",mytpd.intg);
	printf("TPD integral: %.3e\n",mytpd.intg);
	printf("Dose integral: %.3e\n",doseintg);
	printf("Itpd/Idose: %.3f\n",mytpd.intg/doseintg);
  
  #Rate and maximum position
  printf("TPD rate %.2f K/s\n",mytpd.rate);
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
		legend("boxon");
		h=legend(ret.legend);
		set (h, 'fontsize', 10);
	endif;
	
	
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
	plot(dose.t,dose.i*1e10,"linewidth",2,"color",mytpd.color);
	[dosi,endpt]=calculateDoseIntegral(dose,0);
  plot(dose.t(endpt),dose.i(endpt)*1e10,"o")
	printf("Dose integral %.2e\n",dosi);
	
	[maxi,maxidx]=max(dose.i);
	maxt=dose.t(maxidx);
	text(maxt,maxi*1e10,strcat(num2str(mytpd.idx)));
endfunction;


if (index(param.tools,'D'))
	figure(getFigIndex("doseext"));
	hold on;
	ylabel("Dose current (arb.u.)");
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
	ret=iterateTpd(datindex,param,@dispTfitq);
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	
endif;

##################################################
# Log(i) over 1/T plot + linear Eads extraction  #
##################################################
function result=plotInvT(mytpd,param,result);
	mytpd.invT=1./mytpd.T;
	mytpd.logi=log(mytpd.i);
	
	figure(getFigIndex("log"));
	plot(mytpd.invT,mytpd.logi,"linewidth",2,"color",mytpd.color);

	[eads,lnv,win]=findLogEAds(mytpd);
	eVKcalmul=23.0609;
  eVKJmul=96.4869;
	printf("Eads = %f eV, %f kcal/mol, %f kJ/mol\n",eads,eads*eVKcalmul,eads*eVKJmul);
	txtx=mytpd.invT(1);
	txty=mytpd.logi(1);
	text(txtx,txty,strcat("<",num2str(mytpd.idx)));
	
	if(isfield(param,"fig") && isfield(param.fig,"disp"))
		figure(getFigIndex("disp"));
		plot(mytpd.T(win.istart),mytpd.i(win.istart),"o");
		plot(mytpd.T(win.iend),mytpd.i(win.iend),"o");
	endif;
	
endfunction

if (index(param.tools,'l'))
	figure(getFigIndex("log"));
	hold on;
	ylabel("log(i)")
	xlabel("Inverse Temperature (1/T)")
	
	iterateTpd(datindex,param,@plotInvT);
	
	
	
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
	xlabel("Wavelength (cm-1)")
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
  saveFig(param,"disp","desorption");
  saveFig(param,'doses','sticking');
  saveFig(param,'doseext','doses');
  saveFig(param,'log','logplot');
  saveFig(param,'eads','eadsest');
  saveFig(param,'press','pressure');
  saveFig(param,'qipress','qms-pressure');
  saveFig(param,'user','user');
  saveFig(param,'IR','FTIR');
  saveFig(param,'modelediff','prefactor');
  saveFig(param,'Tfitq','Trampqual');
  
  saveFig(param,'fit_1Efit',"VSearch1Efit");
  saveFig(param,'fit_dEfit',"VSearchdEfit");
  saveFig(param,'fit_estv_stds',"VSearchEstd");
  saveFig(param,'fit_finaldE',"VSearchdE");
  
endif;
exit(0);

