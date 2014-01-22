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
		if (isfield(tpd,"iN"))
			cutdat=getMassData(tpd,param.displayT,param.mass);
			cutdat.doseintg=dose.integral;
			cutdat.filename=filename;
			if (length(cutdat.i)>0)
				if (index(param.tools,'n'))
					colors=["black";"cyan";"green";"magenta";"red";"yellow"];
					cutdat.color=colors(idx);
				else
					cutdat.color=getLineColor(cutdat.intg,param.monolayer);
				endif;
				cutdat.idx=++counter;
				if (isfield(tpd,"version"));
					if(tpd.version>=20140120)
						result=feval(funcName,cutdat,param,result,press,dose);
					elseif(tpd.version>=20131025)
						result=feval(funcName,cutdat,param,result,press);
					endif;
				else
					result=feval(funcName,cutdat,param,result);
				endif;
				if (index(param.tools,'C'));
					if (yes_or_no("clear plot?"))
						hold off;
						plot(0,0);
						hold on;
					endif
				endif;
				
			else
				printf("No appropriate data in this file.\n");
			endif
		endif
	end
endfunction