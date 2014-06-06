function result=plotEAds(mytpd,param,result);
	tpd.a=[1e13, 1e15, 1e17, 1e19];
	source("~/octave/constants.m");
	cov=mytpd.intg-cumtrapz(mytpd.t,mytpd.i);
	for ai=1:length(tpd.a);
		E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ ((tpd.a(ai)).*cov));
		plot(cov/param.monolayer,E,"color",mytpd.color);
		if (mytpd.idx==1)
			[maxE,maxEpos]=max(E);
			txt=strcat("<",num2str(mytpd.idx),":",num2str(tpd.a(ai)));
			text(cov(maxEpos),maxE,txt);
		endif;
    ai
    tpd.a(ai)
	end
endfunction