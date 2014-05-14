#Cut all arrays in a structure according to a range of time
function result=cutTime(struct,varargin);
	if (isstruct(varargin{1}))
    timecut=varargin{1};
  elseif (length(varargin)==2)
    timecut.min=varargin{1};
    timecut.max=varargin{2};
  endif;
	numpoints=length(struct.t);
	minidx=max(max(find(struct.t<=timecut.min)),1);
	maxidx=min(numpoints,min(find(struct.t >= timecut.max)));
	
	result=cutIndex(struct,numpoints,minidx,maxidx);
endfunction