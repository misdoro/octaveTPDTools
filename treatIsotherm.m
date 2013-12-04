function result=treatIsotherm(mytpd,param,result);

isot=getIsotherm(mytpd);

figure(1);
hold on;
plot(isot.t-isot.t(1),isot.i,'color',mytpd.color);


figure(2);
hold on;
plot(isot.t-isot.t(1),isot.T,'color',mytpd.color);

endfunction

###################################################
#Get the isothermal part of the desorption curve
####################################################
function result=getIsotherm(mytpd)

#Remove the beginning of the data
	cutT.min=min(mytpd.T)+10;
	cutT.max=max(mytpd.T);
	if (mytpd.T(end)>=cutT.max)
		tpdc=cutTemp(mytpd,cutT);
	else
		tpdc=mytpd;
	endif;

#Find the start of the isotherm
	tpdc.dT=diff(tpdc.T);
	tpdc.dT=[tpdc.dT(1);tpdc.dT];
	tpdc.dTs=supsmu(tpdc.t,tpdc.dT,'span',0.01);

	isi=min(find(tpdc.dTs<0));
#Cut and find the end of the isotherm
	tpdcc=cutIndex(tpdc,length(tpdc.T),isi,length(tpdc.T));
	ise=min(find(tpdcc.dTs>tpdc.dTs(1)))-20;
	if (length(ise)==0)
		result=tpdcc;
	else
		result=cutIndex(tpdcc,length(tpdcc.T),1,ise);
	endif;

endfunction;

