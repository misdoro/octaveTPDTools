#!/usr/bin/octave -qi
#Arguments: mid, min T, max T, display, monolayer

indata.filenames=findDatFiles(pwd);



if (nargin==0)
	printf("\n\
  ################################################################\n\
  # TPD data treatment program, LPMAA/LERMA, UPMC, Paris         #\n\
  ################################################################\n\
  #Usage: display.m act [M1,M2,MN startTemp endTemp monolayer]   #\n\
  #Actions:                                                      #\n\
  #-----------display --------------                             #\n\  
  #d = plot TPDs                                                 #\n\
  #D = plot dose graph                                           #\n\
  #r = plot real non-smooth data                                 #\n\
  #x = extract the baseline from the data                        #\n\
  #c = use coverage-based colors                                 #\n\
  #f = plot the temperature fit quality                          #\n\
  #p = plot pressure                                             #\n\
  #P = plot pressure vs Iqms                                     #\n\
  #R = plot the IR spectrum, add points to the TPD plot          #\n\
  #-----------Ea extract -----------                             #\n\
  #l = log plot of i over 1/T                                    #\n\
  #e = energy estimation using inversion plot over coverage      #\n\
  #t = plot temperature points on the inversion curves every 5K  #\n\
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

param.monolayer=3.2e-09;#Xe /crystal
param.monolayer=1.65e-09;#Xe/HOPG
param.monolayer=2e-08;#Xe /amorph 1,2e-6 A*s 100K
#param.monolayer=1e-5;#H2O
param.monolayer=useArgument(argv(),5,param.monolayer);
param=loadParamFile(param);


#By default, use the first MID and available temperature range from the first input file.
param=getFileInfo(indata,param);
param.displayT.min=useArgument(argv(),3,param.displayT.min);
param.displayT.max=useArgument(argv(),4,param.displayT.max);
param.mass=useArgument(argv(),2,param.mass);
for midx=1:length(param.mass)
	printf("Displayed temperature: %.2f to %.2f K\nDisplayed MID: %d\n"
	,param.displayT.min,param.displayT.max,param.mass(midx));
end

indata.doses=loadDoses(indata.filenames,param.mass(1));

indata.sorted=sortrows(indata.doses,2);


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
		mytpd.i_sm=supsmu(mytpd.T,mytpd.i,'spa',0.005);
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
	legendtext=strcat(txt,":",fn,"(",num2str(mytpd.intg/param.monolayer,"%3.2f"),"ML)");
	
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
  if (index(param.tools,'s'))
    table=[mytpd.t-mytpd.t(1),mytpd.T,mytpd.i,mytpd.pi];
    save("-text",strcat(mytpd.filename,".asc"),"table");
  endif
endfunction

if (index(param.tools,'d'))
	param.fig.disp=++param.figindex;
	figure(param.fig.disp);
	hold on;
	ret=iterateTpd(indata,param,@plotTPD);
	ylabel("Desorption flow (arb.u.)");
	xlabel("Temperature (K)");
	if (isfield(ret,"legend"))
		legend("boxon");
		h=legend(ret.legend);
		set (h, 'fontsize', 10);
	endif;
	
	
	#Plot doses
	if (isfield(ret,"doses"))
		param.fig.doses=++param.figindex;
		figure(param.fig.doses);
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
    param.fig.user=++param.figindex;
    figure(param.fig.user);
	
    #First check if we have the iteration script
    if (init.err>=0);
      source("init.m");
    endif;
    ret=iterateTpd(indata,param,@iterTPD);
    
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
	dosi=calculateDoseIntegral(dose,0);
	printf("Dose integral %.2e\n",dosi);
	
	[maxi,maxidx]=max(dose.i);
	maxt=dose.t(maxidx);
	text(maxt,maxi*1e10,strcat(num2str(mytpd.idx)));
endfunction;


if (index(param.tools,'D'))
	param.fig.doseext=++param.figindex;
	figure(param.fig.doseext);
	hold on;
	ylabel("Dose current (arb.u.)");
	xlabel("Time (s)");
	ret=iterateTpd(indata,param,@plotDoses);
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
	param.fig.press=++param.figindex;
	figure(param.fig.press);
	hold on;
	ret=iterateTpd(indata,param,@plotP);
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
	param.fig.qipress=++param.figindex;
	figure(param.fig.qipress);
	hold off;
	ylabel("Iqms, A");
	xlabel("Pressure (torr)");
	ret=iterateTpd(indata,param,@dispplotqP);
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	
endif;

#################################################
# Plot temperature fit quality
#################################################
if (index(param.tools,'f'))
	param.fig.Tfitq=++param.figindex;
	figure(param.fig.Tfitq);
	hold on;
	ylabel("T ramp fit error, K");
	xlabel("T, K");
	ret=iterateTpd(indata,param,@dispTfitq);
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
	
	figure(param.fig.log);
	plot(mytpd.invT,mytpd.logi,"linewidth",2,"color",mytpd.color);

	[eads,lnv,win]=findLogEAds(mytpd);
	eVKcalmul=23.0609;
  eVKJmul=96.4869;
	printf("Eads = %f eV, %f kcal/mol, %f kJ/mol\n",eads,eads*eVKcalmul,eads*eVKJmul);
	txtx=mytpd.invT(1);
	txty=mytpd.logi(1);
	text(txtx,txty,strcat("<",num2str(mytpd.idx)));
	
	if(isfield(param,"fig") && isfield(param.fig,"disp"))
		figure(param.fig.disp);
		plot(mytpd.T(win.istart),mytpd.i(win.istart),"o");
		plot(mytpd.T(win.iend),mytpd.i(win.iend),"o");
	endif;
	
endfunction

if (index(param.tools,'l'))
	param.fig.log=++param.figindex;
	figure(param.fig.log);
	hold on;
	ylabel("log(i)")
	xlabel("Inverse Temperature (1/T)")
	
	iterateTpd(indata,param,@plotInvT);
	
	
	
endif;


##################################################
#Eads inversion	                          			 #
##################################################
if (index(param.tools,'e'))
	param.fig.eads=++param.figindex;
	figure(param.fig.eads);
	hold on;
	ylabel("Ea (eV)")
	xlabel("Coverage (ML)")
	ret3=iterateTpd(indata,param,@plotEAds);
endif

#######################################################
# TPD modeling, uses (user-defined) function fitModel #
#######################################################

if (index(param.tools,'m'));
  param.fig.model=++param.figindex;
  param.fig.modelfit=++param.figindex;
  param.fig.modelediff=++param.figindex;
	figure(param.fig.model);
	hold on;
	pkg load odepkg;
	iterateTpd(indata,param,@fitModel);
  figure(param.fig.model);
	xlabel("Temperature (K)")
	ylabel("Desorption signal");
  
  #Plot fit results, if any
  fitfiles=findDatFiles(pwd,"fit");
  if (length(fitfiles)>1)
    figure(param.fig.modelediff);
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
	printInfo(indata,param);
endif

########################################################
# Experimental isotherm treatment, not complete        #
########################################################
	
if (index(param.tools,'T'));
	figure(++param.figindex);
	hold on;
	baseparam=iterateTpd(indata,param,@treatIsotherm);
endif

########################################################
# Infra-red spectrum plot, with colour points on TPD   #
########################################################
if (index(param.tools,'R'));
	param.fig.IR=++param.figindex;
  figure(param.fig.IR);
  ylabel("Absorbance")
	xlabel("Wavelength (cm-1)")
	hold on;
  [f.info, f.err, f.msg]=stat("IR.irdat");
	if (f.err>=0);
    load ("IR.irdat");
	  points=iterateTpd(indata,param,@plotIRData,IRDAT);
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
if(isfield(param,"fig"))
	input("Adjust figures if needed, then press enter to save, ctrl-c to exit");
	if (isfield(param.fig,"disp"))
		print(param.fig.disp,"desorption.png","-dpng","-r300");
		printf("Saved desorption image\n");
	endif;
	
	if (isfield(param.fig,"doses"))
		print(param.fig.doses,"sticking.png","-dpng","-r300");
		printf("Saved sticking graph image\n");
	endif;
	if (isfield(param.fig,"doseext"))
		print(param.fig.doseext,"doses.png","-dpng","-r300");
		printf("Saved dose information image\n");
	endif;
	
	if (isfield(param.fig,"log"))
		print(param.fig.log,"logplot.png","-dpng","-r300");
		printf("Saved log fit image\n");
	endif;
	if (isfield(param.fig,"eads"))
		print(param.fig.eads,"eadsest.png","-dpng","-r300");
		printf("Saved Ea inversion image\n");
	endif;
	if (isfield(param.fig,"blf"))
		print(param.fig.blf,"baselinefix.png","-dpng","-r300");
		printf("Saved the baseline-corrected image\n");
	endif;
	if (isfield(param.fig,"blcalib"))
		print(param.fig.blcalib,"baselinecal.png","-dpng","-r300");
		printf("Saved the baseline calibration image\n");
	endif;
	if (isfield(param.fig,"press"))
		print(param.fig.press,"pressure.png","-dpng","-r300");
		printf("Saved the pressure vs T image\n");
	endif;
	if (isfield(param.fig,"qipress"))
		print(param.fig.qipress,"qms-pressure.png","-dpng","-r300");
		printf("Saved the Iqms vs pressure image\n");
	endif;
  if (isfield(param.fig,"user"))
		print(param.fig.user,"user.png","-dpng","-r300");
		printf("Saved the user-specified image\n");
	endif;
  if (isfield(param.fig,"IR"))
		print(param.fig.IR,"FTIR.png","-dpng","-r300");
		printf("Saved the FTIR image\n");
	endif;
  if (isfield(param.fig,"modelediff"))
    print(param.fig.modelediff,"prefactor.png","-dpng","-r300");
		printf("Saved the prefactor fit image\n");
  endif;
  if (isfield(param.fig,"Tfitq"))
    print(param.fig.Tfitq,"Trampqual.png","-dpng","-r300");
		printf("Saved the ramp quality image\n");
  endif;
endif;
exit(0);

