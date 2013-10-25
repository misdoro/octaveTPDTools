function result=calibrateBaseLine(mytpd,param,result,press);
	%Function to calibrate baseline parameters using TPD position dose
	
	result=plotTPD(mytpd,param,result);
	%Initial parameters guess
	par(1)=0;	#delay
	par(2)=min(mytpd.pi); #base
	par(3)=max(mytpd.i)/(max(mytpd.pi)-par(2)); #scale
	
	%Plot initial guess
	pio=interpBaseLine(par,mytpd,press);
	plot(mytpd.T,pio,"color","green",'linewidth',1);
	drawnow();
	
	
	%Fit baseline using fminsearch
	tic();
	fh=@(param)shiftp(param,mytpd,press);
	options(6)=2;
	[paro,mindiff]=fminsearch(fh,par,options,[])
	toc()
	
	%Plot fits and errors
	mytpd.pio=interpBaseLine(paro,mytpd,press);
	tpd_c=cutNA(mytpd,"pio");
	plot(tpd_c.T,tpd_c.pio,"color","blue",'linewidth',1);
	plot(tpd_c.T,(tpd_c.i-tpd_c.pio)*100,"color","black",'linewidth',1);
	drawnow();
	
	result.baselineopt{mytpd.idx}=paro;
endfunction;

function result = shiftp(par,tpd,press)
	tpd.pis=interpBaseLine(par,tpd,press);
	tpd_cut=cutNA(tpd,"pis");
	result=sumsq(tpd_cut.i-tpd_cut.pis);
endfunction;