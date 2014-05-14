#Return a list of .dat files from a given directory
function result=findDatFiles(directory)
	files=readdir(directory);
	len=length(files);

	result={};

	for idx=1:length(files);
		if (regexpi(files{idx},"\\.dat$")>1 && ~isdir(files{idx}))
			result=[result,files{idx}];
		endif
	end;
  