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
				param.selectedmass=param.mass(midx);
				cutdat=getMassData(tpd,param.displayT,param.selectedmass);
			  if (length(cutdat.i)>0)
          if (index(param.tools,'x'))
            cutdat=findBaseLine(cutdat);
            purintg=trapz(cutdat.t,cutdat.ipur);
            intg=trapz(cutdat.t,cutdat.i);
            printf("Baseline integral part: %d percent\n",(1-purintg/intg)*100);
            cutdat.i=cutdat.ipur;
            cutdat.intg=purintg;
          endif 
          cutdat.doseintg=dose.integral;
				  cutdat.filename=filename;
				
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