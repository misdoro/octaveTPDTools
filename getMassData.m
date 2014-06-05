#Extract from the tpd structure data about requested mass in a specified temperature range
function result = getMassData(tpd,Trange,mass);
	massidx=find(tpd.mids==mass);
	zerosi=min(find(tpd.mids>0));
	if (mass>0 && massidx>=1)
		result=tpd;
		colidx=massidx-zerosi+1;
		result.i=tpd.iN(:,colidx);
		if (isstruct(Trange) && Trange.max>0)
			result=findBaseLine(cutTemp(result,Trange));
      purintg=trapz(result.t,result.ipur);
      intg=trapz(result.t,result.i);
      printf("Baseline integral part: %d percent\n",(1-purintg/intg)*100);
      result.i=result.ipur;
		endif;
		result.intg=trapz(result.t,result.i);
		result.mass=mass;
	else
		result.T=[];
		result.i=[];
		result.intg=0;
	endif;
endfunction;