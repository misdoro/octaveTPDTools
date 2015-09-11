#Iterate through sorted TPD files, call function on each one. returns function call result structure.
#Uses:	indata.filenames,
#	indata.ordDoses,
#	param.displayT,
#	param.mass
#Sets:	cutdat.filename
#	cutdat.doseintg
#	cutdat.color
#Calls: function(cutdat,param,prevresult)
function result=iterateTpd(indata,param,funcName,varargin);
	result.idx=[];
  filesCount=rows(indata.ordDoses);
  plotCount=length(param.mass)*filesCount;
	counter=0;
	for i=1:filesCount;
    idx=filesCount+1-i;
		filename=indata.filenames{indata.ordDoses(idx,1)};
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
          if (~isfield(cutdat,"rate"))
            cutdat.rate=0;
          endif
          
          cutdat.baseline=0;
          if (index(param.tools,'x'))
            totintg=trapz(cutdat.t,cutdat.i);
            if (getOptionValue(param,filename,"autobline",0))
              [bl,minidx]=min(cutdat.i);
              imin=max(1,minidx-5);
              imax=min(minidx+5,length(cutdat.i));
              cutdat.baseline=mean(cutdat.i(imin:imax));
              cutdat.i-=cutdat.baseline;
            else
              cutdat.baseline=getOptionValue(param,filename,"bline",0);
              cutdat.i-=cutdat.baseline;
            endif
            cutdat.intg=trapz(cutdat.t,cutdat.i);
            printf("Baseline contribution to TPD integral: %d %%\n",(1-cutdat.intg/totintg)*100);
          endif 
          cutdat.doseintg=dose.integral;
				  cutdat.filename=filename;
					cutdat.idx=++counter;
          cutdat.fileidx=idx;
          
          #Define plot color
					if (index(param.tools,'c'))
            cutdat.color=getLineColor(cutdat.intg,param.monolayer);
					else
            cutdat.color=getLineColor(4*(1-counter/plotCount),1,2);
					endif;
          #check if we have color in param.m for this file
          cutdat.color=getOptionValue(param,filename,"color",cutdat.color);
          
          #Check if we are asked to convert the data in ML/K units
          if (index(param.tools,'N'))
            cutdat.i=cutdat.i./(param.monolayer.*(cutdat.rate));
          endif
          
					result=feval(funcName,cutdat,param,result,press,dose,varargin);
				else
					printf("No appropriate data in file %s for mid %d in T range [%d, %d]\n"
          ,filename,param.mass(midx),param.displayT.min,param.displayT.max);
				endif
			end;
		endif
	end
endfunction