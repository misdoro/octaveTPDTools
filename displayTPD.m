#!/usr/bin/octave --persist
#Arguments: mid, min T, max T, display, monolayer

input.filenames=findDatFiles(pwd);

input.doses=loadDoses(input.filenames);

input.sorted=sortrows(input.doses,2);

if (nargin==0)
	printf("\n\
	#\n\
	#\n\
	#\n\
	#Usage: display.m actions mass startTemp endTemp monolayer \n\
	#Actions:\n\
	#d = plot TPDs\n\
	#s = surf TPDs\n\
	#p = plot pressure\n\
	#c = plot pressure-corrected TPDs\n\
	#C = calibrate pressure correction\n\
	#I = interactive parameters\n\
	#l = log plot of i over 1/T\n\
	#e = energy estimation using inversion plot over coverage\n\
	#m = try to model TPD with 1-st order process\n\
	#i = print files info: available masses and T range\n\
	#\n\
	")
endif

param.monolayer=useArgument(argv(),5,3.2e-09);
#Xe /crystal
#param.monolayer=2e-08;#Xe /amorph 1,2e-6 A*s 100K
#param.monolayer=1e-5;#H2O

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
	plot(mytpd.T,mytpd.i_sm,"linewidth",2,"color",mytpd.color);
	
	[maxi,maxidx]=max(mytpd.i_sm);
	maxT=mytpd.T(maxidx);
	text(maxT,maxi,strcat(num2str(mytpd.idx)));
	
	legendtext=strcat(num2str(mytpd.idx),":",mytpd.filename,"(",num2str(mytpd.intg/param.monolayer,"%3.2f"),"ML)");
	
	result=retAppend(result,"legend",legendtext);
	result=retAppend(result,"doses",mytpd.doseintg);
	result=retAppend(result,"integrals",mytpd.intg);
endfunction

if (index(param.tools,'d'))
	figure(++param.figindex);
	hold on;
	ret=iterateTpd(input,param,@plotTPD);
	ylabel("Desorption flow (arb.u.)");
	xlabel("Temperature (K)");
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	print(param.figindex,"desorption.png","-dpng","-r300");
	
	#Plot doses
	if (isfield(ret,"doses"))
		figure(++param.figindex)
		plotDoses(ret);
		print(param.figindex,"sticking.png","-dpng","-r300");
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
	iterateTpd(input,param,@surfTPD);
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
	ret=iterateTpd(input,param,@plotP);
	ylabel("Pressure, torr");
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
	
	plot(mytpd.invT,mytpd.logi,"linewidth",2,"color",mytpd.color);

	eads=findLogEAds(mytpd);
	eVKcalmul=23.0609;
	printf("Eads = %f eV, %f kcal/mol\n",eads,eads*eVKcalmul);
	
	txtx=mytpd.invT(1);
	txty=mytpd.logi(1);
	text(txtx,txty,strcat("<",num2str(mytpd.idx)));
	
endfunction

if (index(param.tools,'l'))
	figure(++param.figindex);
	hold on;
	ret2=iterateTpd(input,param,@plotInvT);
	ylabel("log(i)")
	xlabel("Inverse Temperature (1/T)")
	print(param.figindex,"logplot.png","-dpng","-r300");
endif;


##################################################
#Stupid Eads inversion				 #
##################################################

function result=plotEAds(mytpd,param,result);
	tpd.a=[1e10, 1e11, 1e12, 1e13, 1e14, 1e15];
	source("~/octave/constants.m");
	cov=mytpd.intg-cumtrapz(mytpd.t,mytpd.i);
	for ai=1:length(tpd.a);
		#E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ ((mytpd.rate*tpd.a(ai)).*cov));
		E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ ((tpd.a(ai)).*cov));
		plot(cov,log(E),"color",mytpd.color);
		text(cov(1),E(1),strcat("<",num2str(mytpd.idx),":",num2str(tpd.a(ai))));
	end
endfunction

if (index(param.tools,'e'))
	figure(++param.figindex);
	hold on;
	ret3=iterateTpd(input,param,@plotEAds);
	ylabel("Eads estimation")
	xlabel("Coverage")
endif

#####################################################
# Basic TPD modeling, should use log fit parameters #
#####################################################

function result=drawModel(mytpd,param,result);
	result=plotTPD(mytpd,param,result);
	isc=mytpd.intg/param.monolayer
	[Tode,theta, p]=modelTPD1(mytpd.T,isc,0.8e14,0.416);
	sig=0.9*p*param.monolayer*mytpd.rate;
	plot(Tode,sig);
	plot(Tode,(mytpd.i-sig),"color","red");
endfunction

if (index(param.tools,'m'));
	figure(++param.figindex);
	hold on;
	pkg load odepkg;
	iterateTpd(input,param,@drawModel);
	xlabel("Temperature (K)")
	ylabel("Desorption signal");
endif

########################################################
# Print file info for all dat files in the current dir #
########################################################

if (index(param.tools,'i'));
	printInfo(input);
endif

########################################################
# Calibrate pressure-qms offset and scaling            #
########################################################
	
if (index(param.tools,'C'));
	figure(++param.figindex);
	hold on;
	baseparam=iterateTpd(input,param,@calibrateBaseLine);
	
endif


drawnow();

