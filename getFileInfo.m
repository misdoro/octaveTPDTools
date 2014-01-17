function param=getFileInfo(input,param);
	filename=input.filenames{input.sorted(1,1)};
	load(filename);
	if (isfield(tpd,"version") && tpd.version>=20131025)
		param.displayT.min=min(tpd.T);
		param.displayT.max=max(tpd.T);
		param.mass=tpd.mids(find(tpd.mids>0))(1);
	else
		param.displayT.min=10;
		param.displayT.max=200;
		param.mass=0;
	endif;
endfunction