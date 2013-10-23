function result = cutIndex(struct,numpoints,minidx,maxidx);
	result.cutidx=[minidx,maxidx];
	
	for [val, key] = struct
		if (length(val)==numpoints)
			newval=val(minidx:maxidx,:);
		else
			newval=val;
		endif
		result=setfield(result,key,newval);
	endfor
	
endfunction;