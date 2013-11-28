#Cut all arrays in a structure according to range of T
function result=cutTemp(struct,cutT);
	numpoints=length(struct.T);
	minidx=max(max(find(struct.T<=cutT.min)),1);
	if (isnan(minidx) || length(minidx)==0)
		minidx=1;
	endif;
	maxidx=min(numpoints,min(find(struct.T >= cutT.max)))
	
	if (isnan(maxidx) || length(maxidx)==0);
		maxidx=numpoints;
	endif;
	
	result=cutIndex(struct,numpoints,minidx,maxidx);
endfunction