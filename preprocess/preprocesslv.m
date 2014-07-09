#!/usr/bin/octave --persist
if (nargin>=2);
	tmppath=argv(){2};
endif;
	
tpd.data=dlmread(strcat(tmppath,"/tpd.txt"));
dose.data=dlmread(strcat(tmppath,"/dose.txt"));
dosep.data=dlmread(strcat(tmppath,"/dosep.txt"));
press.data=dlmread(strcat(tmppath,"/press.txt"));

tfile=fopen(strcat(tmppath,"/time.txt"));
timeline=fgetl(tfile);
fclose(tfile);
tpd.time=strptime(timeline,"#%Y.%m.%d-%H:%M:%S");

tpd.mids=tpd.data(1,:);
tpd.mids=tpd.mids(find(tpd.mids>0));
tpd.data=tpd.data(2:end,:);
dose.mids=dose.data(1,:);
dose.data=dose.data(2:end,:);


tpd.t=tpd.data(:,1);
tpd.i=tpd.data(:,2);
tpd.iN=tpd.data(:,2:end);
press.t=press.data(:,1);
press.T=press.data(:,2);
press.p=press.data(:,3);
dose.t=dose.data(:,1);
dose.i=dose.data(:,2);
dose.iN=dose.data(:,2:end);
dose.T=dosep.data(:,2);

#fit temperature range
fitT.min=min(press.T)+10;
fitT.max=max(press.T)-10;



fit=fitTempSlope(press.t,press.T,fitT);
press.To=press.T;
press.T=polyval(fit.poly,press.t);
press.Tfp=fit.poly;
tpd.rate=fit.poly(1);
tpd.T=polyval(fit.poly,tpd.t);

pkg load optim;

tpd.pi=interp1(press.t,press.p,tpd.t);

#figure(1);
#hold on;
#plot(press.T,press.p,"color","red",'linewidth',2);
#plot(tpd.T,tpd.i,"color","green",'linewidth',2);

#figure(2);
#plot(dose.t,dose.i);

#Detect and remove spikes from TPD and dose
for i=1:size(tpd.iN,2)
  tpd.iN(:,i)=removeSpikes(tpd.t,tpd.iN(:,i));
end


##################
#Dose integrals
for midx=1:length(dose.mids);
	mid=dose.mids(midx);
	[dose.integral(midx),dose.intimax(midx)]=calculateDoseIntegral(getMassData(dose,[],mid));
end;

#Export temperature range
exportT.min=min(press.To)+1;
exportT.max=max(press.To)-1;

tpd_out=cutTemp(tpd,exportT);
tpd_out.integral=trapz(tpd_out.t,tpd_out.i);
for i=1:size(tpd_out.iN,2)
  tpd_out.intg(i)=trapz(tpd_out.t,tpd_out.iN(:,i));
end
tpd=tpd_out;
tpd.model=0;

dose_out.integral=dose.integral;
dose_out.T=mean(dose.T);
dose_out.mids=dose.mids;
dose_out.t=dose.t;
dose_out.iN=dose.iN;
dose=dose_out;




##########################
#Output data truncation and save
if (nargin>=1);
	filename=argv(){1};
	if (!isempty(strfind(filename,".dat" )))
		tpd.version=20140709;
		save("-binary",filename,"tpd","dose","press");
		
		comb=[tpd.t,tpd.T,tpd.i,tpd.pi];
		doseint=dose.integral;
		save("-text",strcat(filename,".asc"),"doseint","comb");
	endif
endif

