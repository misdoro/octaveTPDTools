function retstr=retAppend(retstr,field,value);
%Append value to the array field of the structure retstr.
% Syntax:
% retstr=retAppend(retstr,field,value);
%
	if (isfield(retstr,field))
		oldval=getfield(retstr,field);
		newval=[oldval;value];
		retstr=setfield(retstr,field,newval);
	else
		retstr=setfield(retstr,field,[value]);
	endif
endfunction
