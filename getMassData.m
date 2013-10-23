#Extract from the tpd structure data about requested mass in a specified temperature range
function result = getMassData(tpd,Trange,mass);
	mass
	massidx=find(tpd.mids==mass);
	if (massidx>1)
		result=tpd;
		result.i=tpd.iN(:,massidx-1);
		result=cutTemp(result,Trange);
		result.intg=trapz(result.t,result.i);
	else
		result.T=[];
		result.i=[];
		result.intg=0;
	endif;
endfunction;