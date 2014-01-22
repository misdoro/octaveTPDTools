#Iterate through sorted TPD files, call function on each one. returns function call result structure.
#Uses:	indata.filenames,
#	indata.sorted,
#	param.displayT,
#	param.mass
#Sets:	cutdat.filename
#	cutdat.doseintg
#	cutdat.color
#Calls: function(cutdat,param,prevresult)
function result=iterateTpd(indata,param,funcName);
	result.idx=[];
	counter=0;
	for idx=1:rows(indata.sorted);
		filename=indata.filenames{indata.sorted(idx,1)};
		printf("\n-----------------\n%s\n",filename);
		load(filename);
		if (!exist("dose","var") || !isstruct(dose));
			dose.data=[]
		endif;
		if (!exist("press","var") || !isstruct(press));
			press.t=[]
		endif;
		if (isfield(tpd,"iN"))
			for midx=1:length(param.mass);
				cutdat=getMassData(tpd,param.displayT,param.mass(midx));
				cutdat.doseintg=dose.integral;
				cutdat.filename=filename;
				if (length(cutdat.i)>0)
					cutdat.idx=++counter;
					if (index(param.tools,'n'))
						colors=["black";"cyan";"green";"magenta";"red";"yellow"];
						cutdat.color=colors(counter);
					else
						cutdat.color=getLineColor(cutdat.intg,param.monolayer);
					endif;
					result=feval(funcName,cutdat,param,result,press,dose);
				else
					printf("No appropriate data in file %s for mid %d in T range [%d, %d]\n"\
					,filename,param.mass(midx),param.displayT.min,param.displayT.max);
				endif
			end;
		endif
	end
endfunction