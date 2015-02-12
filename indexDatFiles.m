function result=indexDatFiles(result,param)
%Function to index dat files in directory
result.rates=loadRates(result.filenames);
result.doses=loadDoses(result.filenames,param.mass(1));
result.ordDoses=sortrows(result.doses,2);
result.ordRates=sortrows(result.rates,2);

endfunction

#Load dose integrals from filenames for sorting
function result=loadRates(filenames)
	result=[];
	for idx=1:length(filenames);
		load(filenames{idx});
		result=[result; idx, tpd.rate];
	endfor;
endfunction;

function result=loadDoses(filenames,mass)
	len=length(filenames);
	result=[];
	for idx=1:len;
		load(filenames{idx});
		if (isfield(tpd,"version") && tpd.version>=20140120 && ~tpd.model)
			dosedat=getMassData(dose,[],mass);
			doseint=calculateDoseIntegral(dosedat);
			result=[result; idx, doseint];
		else
			result=[result; idx, dose.integral];
		endif
	end;
endfunction;