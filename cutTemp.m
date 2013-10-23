#Cut all arrays in a structure according to range of T
function result=cutTemp(struct,cutT);
	numpoints=length(struct.T);
	minidx=max(max(find(struct.T<=cutT.min)),1);
	maxidx=min(numpoints,min(find(struct.T >= cutT.max)));
	
	result=cutIndex(struct,numpoints,minidx,maxidx);
endfunction