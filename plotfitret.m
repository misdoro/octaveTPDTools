function ret=plotfitret(filename);
load(filename);
ret.v=fitret(:,1);
ret.Ea=fitret(:,2);
ret.cov=fitret(:,3);
ret.ssq=fitret(:,4);
semilogx(ret.v,ret.Ea,ret.v,ret.ssq/min(ret.ssq));
endfunction;