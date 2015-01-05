function retval=getOptionValue(par,tpdname,option,defval=0)
%Get value from param structure, in order: file-specific, global or default.
  filepref=strsplit(tpdname,"."){1};
  if (isfield(par,filepref))
    parfs=getfield(par,filepref);
    if (isfield(parfs,option))
      retval=getfield(parfs,option);
      return;
    endif
  endif
  if(isfield(par,option))
    retval=getfield(par,option);
    return;
  else
    retval=defval;
  endif;
endfunction;