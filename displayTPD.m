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

#plot TPD
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


	
#Log(i) over 1/T plot
function result=plotInvT(mytpd,param,result);
	plot(1./mytpd.T,log(mytpd.i),"linewidth",2,"color",mytpd.color);
	txtx=1./mytpd.T(1);
	txty=log(mytpd.i(1));
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



function result=drawModel(mytpd,param,result);
	result=plotTPD(mytpd,param,result);
	[Tode,theta, p]=modelTPD1(mytpd.T,mytpd.intg/param.monolayer,2e13,0.156);
	plot(Tode,0.9*p*param.monolayer*mytpd.rate,"color",mytpd.color);
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

