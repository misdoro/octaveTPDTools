#!/usr/bin/octave --persist
if (nargin>=2);
	tmppath=argv(){2};
endif;
	
tpd.data=dlmread(strcat(tmppath,"/tpd.txt"));
press.data=dlmread(strcat(tmppath,"/press.txt"));

tpd.mids=tpd.data(1,:);
tpd.mids=tpd.mids(find(tpd.mids>0));
tpd.data=tpd.data(2:end,:);
#dose.mids=dose.data(1,:);
#dose.data=dose.data(2:end,:);


tpd.t=tpd.data(:,1);
tpd.i=tpd.data(:,2);
tpd.iN=tpd.data(:,2:end);
press.t=press.data(:,1);
press.T=press.data(:,2);
press.p=press.data(:,3);
#dose.t=dose.data(:,1);
#dose.i=dose.data(:,2);
#dose.iN=dose.data(:,2:end);
#dose.T=dosep.data(:,2);

#fit temperature range
#fitT.min=min(press.T)+10;
#fitT.max=max(press.T)-10;



#fit=fitTempSlope(press.t,press.T,fitT);
press.To=press.T;
#press.T=fit.points;
#press.Tfp=fit.poly;
#tpd.rate=fit.poly(1);
#tpd.T=polyval(fit.poly,tpd.t);
tpd.T=interp1(press.t,press.T,tpd.t);

pkg load optim;

tpd.pi=interp1(press.t,press.p,tpd.t);

#figure(1);
#hold on;
#plot(press.T,press.p,"color","red",'linewidth',2);
#plot(tpd.T,tpd.i,"color","green",'linewidth',2);

#figure(2);
#plot(dose.t,dose.i);


##################
#Dose integrals
#for midx=1:length(dose.mids);
#	mid=dose.mids(midx);
#	[dose.integral(midx),dose.intimax(midx)]=calculateDoseIntegral(getMassData(dose,[],mid));
#end;

#Export temperature range
#exportT.min=min(press.To)+1;
#exportT.max=max(press.To)-1;

#tpd_out=cutTemp(tpd,exportT);
tpd.integral=trapz(tpd.t,tpd.i);
#tpd=tpd_out;
tpd.idx=1;
tpd.color=[0,0,0];
tpd.model=0;

#dose_out.integral=dose.integral;
#dose_out.T=mean(dose.T);
#dose_out.mids=dose.mids;
#dose_out.t=dose.t;
#dose_out.iN=dose.iN;
#dose=dose_out;
dose.T=[];
dose.mids=[];
dose.integral=0;

##########################
#Output data truncation and save
if (nargin>=1);
	filename=argv(){1};
	if (!isempty(strfind(filename,".dat" )))
		tpd.version=20140207;
		save("-binary",filename,"tpd","dose","press");
		
		comb=[tpd.t,tpd.T,tpd.i,tpd.pi];
		doseint=0;
		save("-text",strcat(filename,".asc"),"doseint","comb");
	endif
endif

