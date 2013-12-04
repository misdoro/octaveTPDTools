#Cut all arrays in a structure according to range of T
function result=cutTemp(struct,cutT);
	numpoints=length(struct.T);
	minidx=max(max(find(struct.T<=cutT.min)),1);
	if (isnan(minidx) || length(minidx)==0)
		minidx=1;
	endif;
	
	maxidx=min(numpoints,min(find(struct.T >= cutT.max)));
	
	if (isnan(maxidx) || length(maxidx)==0);
		maxidx=numpoints;
	endif;
	
	if (minidx>maxidx)
		minidx=max(max(find(struct.T(1:maxidx)<=cutT.min)),1);
	endif;
	
	if (minidx==numpoints || maxidx==1)
		result.T=[];
		result.t=[];
		result.i=[];
		result.intg=0;
		return;
	endif;
	
	result=cutIndex(struct,numpoints,minidx,maxidx);
endfunction