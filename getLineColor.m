#Get the line color as a function of coverage integral
function retval=getLineColor(integral,monolay=1)
	if (integral<=0)
		integral=0;
	endif;
	imrel=integral/monolay;
	if (integral>monolay)
		r=1;
		g=min(1,log(integral/monolay));
	else
		r=(imrel);
		g=0;
	endif;
	b=1-r;
	retval=[r,g,b];
endfunction