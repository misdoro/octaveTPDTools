function ret=invertEa(mytpd,vi)
  constants;
  if (~isfield(mytpd,"cov"))
    mytpd.cov=trapz(mytpd.t,mytpd.i)-cumtrapz(mytpd.t,mytpd.i);
  endif
  [maxd,maxdi]=max(mytpd.i);
  
  ret.E=-(R_eV/Na).* mytpd.T .*log(mytpd.i ./ (vi.*(mytpd.cov+eps)));
    
  if(isinf(ret.E(1)))
    ret.startidx=max(find(isinf(ret.E(1:maxdi))));
  else
    ret.startidx=20;
  endif
  ret.endidx=min([find(isnan(ret.E));length(ret.E)-10]);
endfunction