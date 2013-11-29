function result=treatIsotherm(mytpd,param,result);

#Remove the beginning of the curve
cutT.min=min(mytpd.T)+10;
cutT.max=max(mytpd.T);
tpdc=cutTemp(mytpd,cutT);

#Find the start of the isotherm
tpdc.dT=diff(tpdc.T);
tpdc.dT=[tpdc.dT(1);tpdc.dT];
tpdc.dTs=supsmu(tpdc.t,tpdc.dT,'span',0.01);

isi=min(find(tpdc.dTs<0));

tpdcc=cutIndex(tpdc,length(tpdc.T),isi,length(tpdc.T));
ise=min(find(tpdcc.dTs>tpdc.dTs(1)))-20;

isot=cutIndex(tpdcc,length(tpdcc.T),1,ise);
plot(isot.t-isot.t(1),isot.i);


endfunction