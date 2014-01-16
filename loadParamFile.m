function ret=loadParamFile(par);
	[f.info, f.err, f.msg]=stat("param.m");
	if (f.err>=0);
		printf("Loading param\n");
		source("param.m");
	endif;
	ret=par;
endfunction;