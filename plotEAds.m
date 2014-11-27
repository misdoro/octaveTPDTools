function result=plotEAds(mytpd,param,result);
  par.vs=[1e15,1e17,1e19];
  par=loadParamFile(par);
	source("~/octave/constants.m");
	mytpd.cov=trapz(mytpd.t,mytpd.i)-cumtrapz(mytpd.t,mytpd.i);
  cml=mytpd.cov/param.monolayer;
  
	for vi=par.vs;
		inv=invertEa(mytpd,vi);
    E=inv.E;
    maxbeg=inv.startidx;
    mini=inv.endidx;
    plot(cml(maxbeg:mini),E(maxbeg:mini),"color",mytpd.color);
		if (mytpd.idx==1)
			txt=sprintf("< %.1e ",vi);
			text(cml(mini-1),E(mini-1),txt);
		endif;
    sprintf("paramtools:%s",param.tools);
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