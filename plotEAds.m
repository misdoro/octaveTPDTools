function result=plotEAds(mytpd,param,result);
  par.vs=[1e15,1e17,1e19];
  par=loadParamFile(par);
	tpd.a=par.vs;
	source("~/octave/constants.m");
	cov=trapz(mytpd.t,mytpd.i)-cumtrapz(mytpd.t,mytpd.i);
  [maxd,maxdi]=max(mytpd.i);
	for ai=1:length(tpd.a);
		E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ ((tpd.a(ai)).*(cov+eps)));
    cml=cov/param.monolayer;
    if(isinf(E(1)))
      maxbeg=max(find(isinf(E(1:maxdi))))
    else
      maxbeg=20;
    endif
    
    mini=min([find(isnan(E));length(E)-10]);
		plot(cml(maxbeg:mini),E(maxbeg:mini),"color",mytpd.color);
		if (mytpd.idx==1)
			txt=sprintf("< %.1e ",tpd.a(ai));
			text(cml(mini-1),E(mini-1),txt);
		endif;
    if (index(param.tools,'t'))
      minT=5*ceil(min(mytpd.T)/5);
      maxT=5*floor(max(mytpd.T)/5);
      numstep=1+(maxT-minT)/5;
      temps=linspace(minT,maxT,numstep);
      for i=1:numstep
        tidx=max(find(mytpd.T<=temps(i)));
        if (tidx>maxbeg && tidx<mini)
          txt=sprintf("< %dK",temps(i));
          text(cml(tidx),E(tidx),txt);
        endif
      endfor;
    endif;
	end
endfunction