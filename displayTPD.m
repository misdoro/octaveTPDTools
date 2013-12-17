#!/usr/bin/octave --q
#Arguments: mid, min T, max T, display, monolayer

indata.filenames=findDatFiles(pwd);

indata.doses=loadDoses(indata.filenames);

indata.sorted=sortrows(indata.doses,2);

if (nargin==0)
	printf("\n\
	#\n\
	#\n\
	#\n\
	#Usage: display.m actions mass startTemp endTemp monolayer \n\
	#Actions:\n\
	#d = plot TPDs\n\
	#n = use index-based colors\n\
	#s = surf TPDs\n\
	#p = plot pressure\n\
	#c = plot pressure-corrected TPDs\n\
	#C = calibrate pressure correction\n\
	#I = interactive parameters\n\
	#l = log plot of i over 1/T\n\
	#e = energy estimation using inversion plot over coverage\n\
	#m = try to model TPD with 1-st order process\n\
	#i = print files info: available masses and T range\n\
	#v = ask to clear plot after each iteration\n\
	#T = treat isotherm desorption\n\
	#\n\
	")
endif

param.monolayer=3.2e-09;#Xe /crystal
param.monolayer=1.65e-09;#Xe/HOPG
param.monolayer=2e-08;#Xe /amorph 1,2e-6 A*s 100K
#param.monolayer=1e-5;#H2O
param.monolayer=useArgument(argv(),5,param.monolayer);
param.displayT.min=useArgument(argv(),3,45);
param.displayT.max=useArgument(argv(),4,70);
param.mass=useArgument(argv(),2,130);

param.tools=useArgument(argv(),1,"i");
param.figindex=0;%Current figure index, to be autoincremented by processings

pkg load optim;

#################################################
# Plot TPD
#################################################
function result=plotTPD(mytpd,param,result);
	
	mytpd.i_sm=supsmu(mytpd.T,mytpd.i,'spa',0.005);
	ls="-";
	if (isfield(mytpd,"model"));
		ls="--";
	endif;
	plot(mytpd.T,mytpd.i_sm,"linewidth",2,"linestyle",ls,"color",mytpd.color);
	
	[maxi,maxidx]=max(mytpd.i_sm);
	maxT=mytpd.T(maxidx);
	text(maxT,maxi,strcat(num2str(mytpd.idx)));
	
	legendtext=strcat(num2str(mytpd.idx),":",mytpd.filename,"(",num2str(mytpd.intg/param.monolayer,"%3.2f"),"ML)");
	
	result=retAppend(result,"legend",legendtext);
	result=retAppend(result,"doses",mytpd.doseintg);
	result=retAppend(result,"integrals",mytpd.intg);
	printf("TPD integral: %.3e\n",mytpd.intg);
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
		legend(ret.legend);
	endif;
	
	
	#Plot doses
	if (isfield(ret,"doses"))
		param.fig.doses=++param.figindex;
		figure(param.fig.doses);
		plotDoses(ret);
		print(param.fig.doses,"sticking.png","-dpng","-r300");
	endif;
endif;


#################################################
# Surf TPD
#################################################
function result=surfTPD(mytpd,param,result);
	length(mytpd.mids(2:end))
	length(mytpd.T)
	size(mytpd.iN)
	surfc(mytpd.mids(2:end),mytpd.T,mytpd.iN,'Edgecolor', 'none');
endfunction

if (index(param.tools,'s'));
	figure(++param.figindex);
	hold on;
	iterateTpd(indata,param,@surfTPD);
	ylabel("Desorption flow (arb.u.)");
	xlabel("Temperature (K)");
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
	figure(++param.figindex);
	hold on;
	ret=iterateTpd(indata,param,@plotP);
	ylabel("Pressure (torr)");
	xlabel("Temperature (K)");
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	print(param.figindex,"pressure.png","-dpng","-r300");
	
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
	tpd.a=[1e11, 1e13, 1e15];
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

#####################################################
# Basic TPD modeling, should use log fit parameters #
#####################################################

function result=drawModel(mytpd,param,result);
	result=plotTPD(mytpd,param,result);
	isc=mytpd.intg/param.monolayer
	#[Tode,theta, p]=modelTPD1(mytpd.T,isc,0.8e14,0.416);#H2O params
	#sig=0.9*p*param.monolayer*mytpd.rate
	[Tode,theta, p]=modelTPD1(mytpd.T,5,2e14,0.174,mytpd.rate);#H2O params
	sig=p*3e-09;
	#mytpdd=decimateTPD(mytpd,10);
	#s3=max([0,isc-1])
	#s3=2.5
	#s2=min([0.4,isc*0.4])
	#s1=min([0.6,isc*0.6])
	#[Tode,theta, p1]=modelTPD1(mytpd.T,s1,2e13,0.201,mytpd.rate);#Xe params
	#[Tode,theta, p2]=modelTPD1(mytpd.T,s2,2e13,0.187,mytpd.rate);#Xe params
	#[Tode,theta, p3]=modelTPD1(mytpd.T,s3,1e14,0.14,mytpd.rate);#Xe params
	#sig=(p1*s1+p2*s2+p3)*param.monolayer*mytpd.rate;
	plot(Tode,sig);
	plot(Tode,(mytpd.i-sig),"color","red");
endfunction

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
function result=extractBaseLine(mytpd,param,result,press);
	result=plotTPD(mytpd,param,result);
	mytpd.baseLine=interpBaseLine(param.baseLine,mytpd,press);
	mytpd.i=mytpd.i-mytpd.baseLine;
	printf("Corrected TPD integral %.3e\n",trapz(mytpd.t, mytpd.i));
	mytpd.color='black';
	result=plotTPD(mytpd,param,result);
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
endif;
exit(0);

