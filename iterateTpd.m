#Iterate through sorted TPD files, call function on each one. returns function call result structure.
#Uses:	input.filenames,
#	input.sorted,
#	param.displayT,
#	param.mass
#Sets:	cutdat.filename
#	cutdat.doseintg
#	cutdat.color
#Calls: function(cutdat,param,prevresult)
function result=iterateTpd(input,param,funcName);
	result.idx=[];
	counter=0;
	for idx=1:rows(input.sorted);
		filename=input.filenames{input.sorted(idx,1)};
		printf("\n-----------------\n%s\n",filename);
		load(filename);
		if (isfield(tpd,"iN"))
			cutdat=getMassData(tpd,param.displayT,param.mass);
			cutdat.doseintg=dose.integral;
			cutdat.filename=filename;
			if (length(cutdat.i)>0)
				cutdat.color=getLineColor(cutdat.intg,param.monolayer);
				cutdat.idx=++counter;
				result=feval(funcName,cutdat,param,result);
			else
				printf("No appropriate data in this file.\n");
			endif
		endif
	end
endfunction