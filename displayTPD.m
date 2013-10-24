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
	#Usage: display.m mass startTemp endTemp processings monolayer \n\
	#Processings:\n\
	#d = plot TPDs\n\
	#p = plot pressure\n\
	#c = plot pressure-corrected TPDs\n\
	#C = calibrate pressure correction\n\
	#l = log plot of i over 1/T\n\
	#e = energy estimation using inversion plot over coverage\n\
	#m = try to model TPD with 1-st order process\n\
	#\n\
	#\n\
	")
endif

param.monolayer=useArgument(argv(),5,3.2e-09);
#Xe /crystal
#param.monolayer=2e-08;#Xe /amorph 1,2e-6 A*s 100K
#param.monolayer=1e-5;#H2O

param.displayT.min=useArgument(argv(),2,45);
param.displayT.max=useArgument(argv(),3,70);
param.mass=useArgument(argv(),1,130);

param.tools=useArgument(argv(),4,"d");

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
	figure(1);
	hold on;
	ret=iterateTpd(input,param,@plotTPD);
	ylabel("Desorption flow (arb.u.)");
	xlabel("Temperature (K)");
	if (isfield(ret,"legend"))
		legend("boxon");
		legend(ret.legend);
	endif;
	print(1,"desorption.png","-dpng","-r300");
	
	#Plot doses
	if (isfield(ret,"doses"))
		figure(2)
		plotDoses(ret);
		print(2,"sticking.png","-dpng","-r300");
	endif;
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
	figure(3);
	hold on;
	ret2=iterateTpd(input,param,@plotInvT);
	ylabel("log(i)")
	xlabel("Inverse Temperature (1/T)")
	print(3,"logplot.png","-dpng","-r300");
endif;


##################################################
#Stupid Eads inversion				 #
##################################################

function result=plotEAds(mytpd,param,result);
	#tpd.a=[1e10, 1e11, 1e12, 1e13, 1e14, 1e15];
	tpd.a=1e13;
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
	figure(4);
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
	plot(Tode,sig,"color",mytpd.color);
	plot(Tode,(mytpd.i-sig));
endfunction

if (index(param.tools,'m'));
figure(5);
hold on;
pkg load odepkg;
iterateTpd(input,param,@drawModel);
xlabel("Temperature (K)")
ylabel("Desorption signal");
endif

drawnow();

