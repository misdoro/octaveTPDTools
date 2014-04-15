#Extract from the tpd structure data about requested mass in a specified temperature range
function result = getMassData(tpd,Trange,mass);
	massidx=find(tpd.mids==mass);
	zerosi=min(find(tpd.mids>0));
	if (mass>0 && massidx>=1)
		result=tpd;
		colidx=massidx-zerosi+1;
		result.i=tpd.iN(:,colidx);
		if (isstruct(Trange) && Trange.max>0)
			result=cutTemp(result,Trange);
		endif;
		result.intg=trapz(result.t,result.i);
		result.mass=mass;
	else
		result.T=[];
		result.i=[];
		result.intg=0;
	endif;
endfunction;