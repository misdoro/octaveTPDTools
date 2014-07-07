#Get the line color as a function of coverage integral
function retval=getLineColor(integral,monolay=1,mode=1)
	if (integral<=0)
		integral=0;
	endif;
	imrel=integral/monolay;
  if (mode==1)
	  if (integral>monolay)
  		r=1;
  		g=min(1,log(imrel));
    else
  		r=(imrel);
    	g=0;
    endif;
  	b=1-r;
  else
    if (imrel>3)
      r=1;
      g=min(1,log(imrel-2));
      b=1-r;
    elseif (imrel>2)
      r=imrel-2;
      g=0;
      b=1-r;
    elseif(imrel>1)
      r=0;
      b=imrel-1;
      g=1-b;
    else
      r=0;
      b=0;
      g=imrel;
    endif
  endif
	retval=[r,g,b];
endfunction