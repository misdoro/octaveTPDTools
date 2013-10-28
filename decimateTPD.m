function result=decimateTPD(struct,factor);
%Decimate the incoming structure
	result.decimatefactor=factor;
	pkg load signal;
	numpoints=length(struct.t);
	for [val,key] = struct
		if (rows(val)==numpoints)
			newval=downsample(val,factor);
		else
			newval=val;
		endif
		result=setfield(result,key,newval);
	endfor
endfunction;