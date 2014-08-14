#Return a list of .dat files from a given directory
function result=findDatFiles(directory, extension="dat")
	files=readdir(directory);
	len=length(files);
  regex=sprintf("\\.%s$",extension);
	result={};

	for idx=1:length(files);
		if (regexpi(files{idx},regex)>1 && ~isdir(files{idx}))
			result=[result,files{idx}];
		endif
	end;
  