#Load dose integrals from filenames for sorting
function result=loadDoses(filenames)
	len=length(filenames);
	result=[];
	for idx=1:len;
		load(filenames{idx});
		result=[result; idx, dose.integral];
	end;