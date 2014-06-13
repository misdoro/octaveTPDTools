function result=plotEAds(mytpd,param,result);
	tpd.a=[1e11,1e13, 1e15, 1e17];
	source("~/octave/constants.m");
	cov=trapz(mytpd.t,mytpd.i)-cumtrapz(mytpd.t,mytpd.i);
	for ai=1:length(tpd.a);
		E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ ((tpd.a(ai)).*(cov+eps)));
    cml=cov/param.monolayer;
    if(isinf(E(1)))
      maxbeg=max(find(isinf(E)))
    else
      maxbeg=20;
    endif
    
    mini=min([find(isnan(E));length(E)-10]);
    mini
		plot(cml(maxbeg:mini),E(maxbeg:mini),"color",mytpd.color);
		if (mytpd.idx==1)
			
			txt=sprintf("< %.1e ",tpd.a(ai));
			text(cov(mini-1),E(mini-1),txt);
		endif;
	end
endfunction