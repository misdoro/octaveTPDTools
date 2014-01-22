#Extract from the tpd structure data about requested mass in a specified temperature range
function result = getMassData(tpd,Trange,mass);
	massidx=find(tpd.mids==mass);
	zerosi=min(find(tpd.mids>0));
	if (massidx>=1)
		result=tpd;
		result.i=tpd.iN(:,massidx-zerosi+1);
		if (isstruct(Trange))
			result=cutTemp(result,Trange);
		endif;
		result.intg=trapz(result.t,result.i);
	else
		result.T=[];
		result.i=[];
		result.intg=0;
	endif;
endfunction;