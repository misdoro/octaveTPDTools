function result=calibrateBaseLine(mytpd,param,result,press);
	%Function to calibrate baseline parameters using TPD position dose
	
	result=plotTPD(mytpd,param,result);
	
	if (index(param.tools,'I'));% Interactive parameters
		par(1) = input("Input initial delay\n");
		par(2) = input("Input pressure base\n");
		par(3) = input("Input pressure scale\n");
	else
		%Initial parameters guess
		[maxi,maxiidx]=max(mytpd.i);
		[maxp,maxpidx]=max(mytpd.pi);
		minp=min(mytpd.pi);
		mini=min(mytpd.i);
		par(1)=mytpd.t(maxpidx)-mytpd.t(maxiidx);
		par(3)=maxi/maxp; #scale
		par(2)=minp-mini/par(3); #base
		%Iterate a bit to make the guess better
		for i=1:3
			par(3)=mean([par(3),maxi/(maxp-par(2))]); #scale
			par(2)=mean([par(2),minp-mini/par(3)]); #base
		endfor;
	endif;
	
	par
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
	
	if (yes_or_no("clear plot?"))
		hold off;
		plot(0,0);
		hold on;
	endif;
	
	result.baselineopt{mytpd.idx}=paro;
endfunction;

function result = shiftp(par,tpd,press)
	tpd.pis=interpBaseLine(par,tpd,press);
	tpd_cut=cutNA(tpd,"pis");
	result=sumsq(tpd_cut.i-tpd_cut.pis);
endfunction;