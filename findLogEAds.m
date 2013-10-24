function eads=findLogEAds(mytpd)
%
% Find Eads from linear regression of log(i) over 1/T plot
% Algorithm may well be improved!
	%Start up to max(i) point
	[maxd,maxdi]=max(mytpd.i);
	
	win.istart=1;
	win.iend=maxdi;
	win.np=maxdi;
	win.bstep=10;
	win.minw=200;
	delta=0.005;
	rmsv=1;
	spvec=[];
	% Find the max error point, use it to cut the fit window
	while(win.np>win.minw && rmsv>delta)
		[sp,rmsv,pf]=findSplitPoint(mytpd,win,0);
		win=stepWin(win,sp);
		spvec=[spvec,sp];
	end
	
	[sp,rmsv]=findSplitPoint(mytpd,win,1);
	plot(mytpd.invT(spvec),mytpd.logi(spvec),"x");
	drawnow();
	
	R_eV=5.189e19;
	Na=6.022141e23;
	eVKcalmul=23.0609;
	eads=(-R_eV*pf(1))/Na;
	
endfunction;

function [maxei,rms,pf]=findSplitPoint(mytpd,win,doplot=false)
%Find a split point for linear fit

	ctpd=cutIndex(mytpd,length(mytpd.T),win.istart,win.iend);
	pf=polyfit(ctpd.invT,ctpd.logi,1);
	pv=polyval(pf,ctpd.invT);
	sqe=(pv-ctpd.logi).^2;
	[maxe,maxei]=max(sqe);
	if(doplot>0);
		plot(ctpd.invT,pv,"color","black",ctpd.invT,sqe*10,"color","red");
	endif
	maxei+=win.istart-1;
	rms=sum(sqe);
endfunction;

function win=stepWin(win,sp);
	istart=win.istart;
	iend=win.iend;
	bstep=win.bstep;
	if (sp<=mean([istart,iend]))
		if (sp<=istart+bstep)
			istart+=bstep;
		else
			istart=sp;
		endif
	else
		if (sp>=(iend-bstep));
			iend-=bstep;
		else
			iend=sp;
		endif;
	endif;
	win.istart=istart;
	win.iend=iend;
	win.np=iend-istart;
endfunction;