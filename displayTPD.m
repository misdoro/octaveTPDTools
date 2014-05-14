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
	#d = plot TPDs                                                 #\n\
	#D = plot dose metrics                                         #\n\
	#n = use index-based colors                                    #\n\
	#p = plot pressure                                             #\n\
	#P = plot pressure vs Iqms                                     #\n\
	#c = plot pressure-corrected TPDs                              #\n\
	#C = calibrate pressure correction                             #\n\
	#I = interactive parameters                                    #\n\
	#l = log plot of i over 1/T                                    #\n\
	#e = energy estimation using inversion plot over coverage      #\n\
	#m = try to model TPD with 1-st order process                  #\n\
	#i = print files info: available masses and T range            #\n\
	#v = ask to clear plot after each iteration                    #\n\
	#T = treat isotherm desorption                                 #\n\
	#x = extract the min(i) from the data                          #\n\
	################################################################\n\n");
endif

param.monolayer=3.2e-09;#Xe /crystal
param.monolayer=1.65e-09;#Xe/HOPG
param.monolayer=2e-08;#Xe /amorph 1,2e-6 A*s 100K
#param.monolayer=1e-5;#H2O
param.monolayer=useArgument(argv(),5,param.monolayer);



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
	mytpd.i_sm=supsmu(mytpd.T,mytpd.i,'spa',0.005);
	
	ls="-";
	if (isfield(mytpd,"model")&& mytpd.model>0);
		ls=":";
	endif;
	if (index(param.tools,'x'))
		mytpd.i_sm-=mini;
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
	if (isfield(mytpd,"version") && mytpd.version>=20140120)
		doseintg=calculateDoseIntegral(getMassData(dose,[],param.selectedmass));
	else
		doseintg=mytpd.doseintg;
	endif
	result=retAppend(result,"doses",doseintg);
	result=retAppend(result,"integrals",mytpd.intg);
	printf("TPD integral: %.3e\n",mytpd.intg);
	printf("Dose integral: %.3e\n",doseintg);
	printf("Itpd/Idose: %.3f\n",mytpd.intg/doseintg);
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
		set (h, 'fontsize', 12);
	endif;
	
	
	#Plot doses
	if (isfield(ret,"doses"))
		param.fig.doses=++param.figindex;
		figure(param.fig.doses);
		plotDoses(ret);
	endif;
endif;

#################################################
# Plot extended dose information
#################################################
function result=plotDoses(mytpd,param,result,press,dose);
	dose=getMassData(dose,[],param.selectedmass);
	plot(dose.t,dose.i,"linewidth",2,"color",mytpd.color)
	[maxi,maxidx]=max(dose.i);
	maxt=dose.t(maxidx);
	text(maxt,maxi,strcat(num2str(mytpd.idx)));
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
	printf("Eads = %f eV, %f kcal/mol\n",eads,eads*eVKcalmul);
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
#Stupid Eads inversion				 #
##################################################

function result=plotEAds(mytpd,param,result);
	tpd.a=[1e13, 1e15, 1e17, 1e19];
	source("~/octave/constants.m");
	cov=mytpd.intg-cumtrapz(mytpd.t,mytpd.i);
	for ai=1:length(tpd.a);
		E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ ((tpd.a(ai)).*cov));
		plot(cov/param.monolayer,E,"color",mytpd.color);
		if (mytpd.idx==1)
			[maxE,maxEpos]=max(E);
			txt=strcat("<",num2str(mytpd.idx),":",num2str(tpd.a(ai)));
			text(cov(maxEpos),maxE,txt);
		endif;
	end
endfunction

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
	figure(++param.figindex);
	hold on;
	pkg load odepkg;
	iterateTpd(indata,param,@fitModel);
	xlabel("Temperature (K)")
	ylabel("Desorption signal");
endif

########################################################
# Print file info for all dat files in the current dir #
########################################################

if (index(param.tools,'i'));
	printInfo(indata,param);
endif

########################################################
# Calibrate pressure-qms offset and scaling            #
########################################################
	
if (index(param.tools,'C'));
	param.fig.blcalib=++param.figindex;
	figure(param.fig.blcalib);
	hold on;
	baseparam=iterateTpd(indata,param,@calibrateBaseLine);
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

if (index(param.tools,'c'));
	param.fig.blf=++param.figindex;
	figure(param.fig.blf);
	hold on;
	param.baseLine=[0,0,0];
	if (index(param.tools,'I'));% Interactive parameters
		printf("Please provide the baseline parameters. \n\
To obtain them, do the TPD in position dose and run displayTPD C mass startT endT\n");
		
		param.baseLine(1) = input("Input time delay");
		param.baseLine(2) = input("Input pressure base\n");
		param.baseLine(3) = input("Input pressure scale\n");
	else
		#baselines(MID,1:3)=[dt,p0,p1]
		#baselines(130,1:3)=[2.1123e+00,   1.8871e-10,   2.0016e-02];
		baselines(130,1:3)=[2.9660e+00   1.1547e-10   2.0236e-02];
		baselines(84,1:3)=[-1.1146e+01   1.0841e-10   8.2433e-02];
		baselines(84,1:3)=[2.5799e+00   7.2170e-11   5.8502e-02];
		baselines(40,1:3)=[-6.4   2.6e-10   2.3516e-01];
		param.baseLine=baselines(param.mass,:);
	endif;
	
	fixedbase=iterateTpd(indata,param,@extractBaseLine);
	
	
endif

drawnow();


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
		print(param.fig.press,"qms-pressure.png","-dpng","-r300");
		printf("Saved the Iqms vs pressure image\n");
	endif;

endif;
exit(0);

