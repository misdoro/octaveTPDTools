function ret=loadFitFile();
  [f.info, f.err, f.msg]=stat("fit.par");
  if (f.err>=0);
    load "fit.par";
  endif;
  if (exist("fitpar","var"))
    ret=fitpar;
  else
    ret=struct();
  endif
endfunction;