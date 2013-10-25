function ret=interpBaseLine(par,tpd,press)
	ret=(interp1(press.t,press.p,tpd.t+par(1)).- par(2) ).* par(3);
endfunction;