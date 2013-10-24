function eads=findLogEAds(mytpd)
%
% Find Eads from linear regression of log(i) over 1/T plot
% Algorithm may well be improved!
	%Start up to max(i) point
	[maxd,maxdi]=max(mytpd.i);
	
	win.istart=1;
	win.iend=maxdi;
	win.iendmax=maxdi;
	win.np=maxdi;
	win.bstep=10;
	win.minw=200;
	delta=0.0001;
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
	[pf,S]=polyfit(ctpd.invT,ctpd.logi,1);
	sqe=(S.yf-ctpd.logi).^2;
	rms=S.normr^2;
	[maxe,maxei]=max(sqe);
	if(doplot>0);
		plot(ctpd.invT,S.yf,"color","black",ctpd.invT,sqe*1000,"color","red");
	endif
	maxei+=win.istart-1;
endfunction;

function win=stepWin(win,sp);
%Adjust window size using split point
	istart=win.istart;
	iend=win.iend;
	bstep=win.bstep;
	if (sp<=mean([istart,iend]))
		if (sp<=istart+bstep)%Max error on the edge of the window, move it!
			istart+=bstep;
			iend=min(iend+bstep-1,win.iendmax);
		else
			istart=sp;
		endif
	else
		if (sp>=(iend-bstep));%Max error on the edge of the window, move it!
			iend-=bstep;
			istart=max(1,istart-bstep+1);
		else
			iend=sp;
		endif;
	endif;
	win.istart=istart;
	win.iend=iend;
	win.np=iend-istart;
endfunction;