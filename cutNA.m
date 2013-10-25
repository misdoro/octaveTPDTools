#Cut all arrays in a structure according to range of time
function result=cutNA(struct,varname);
	cutvar = getfield(struct,varname);
	nans=isnan(cutvar);
	minidx=min(find(nans==0));
	maxidx=max(find(nans==0));
	numpoints=length(struct.t);
	
	result=cutIndex(struct,numpoints,minidx,maxidx);
endfunction