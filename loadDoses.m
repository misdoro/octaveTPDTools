#Load dose integrals from filenames for sorting
function result=loadDoses(filenames,mass)
	len=length(filenames);
	result=[];
	for idx=1:len;
		load(filenames{idx});
		if (isfield(tpd,"version") && tpd.version>=20140120)
			dosedat=getMassData(dose,[],mass);
			doseint=calculateDoseIntegral(dosedat);
			result=[result; idx, doseint];
		else
			result=[result; idx, dose.integral];
		endif
	end;