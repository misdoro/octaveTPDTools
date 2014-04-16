#Cut all arrays in a structure according to a range of time
function result=cutTime(struct,timecut);
	
	numpoints=length(struct.t);
	minidx=max(max(find(struct.t<=timecut.min)),1);
	maxidx=min(numpoints,find(struct.t >= timecut.max));
	
	result=cutIndex(struct,numpoints,minidx,maxidx);
endfunction