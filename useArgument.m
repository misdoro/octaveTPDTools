function ret=useArgument(argv,index,defval)
%
%return numeric value of an argument from argv at specified index, or default value if no such argument is provided
%SYNTAX: ret=useArgument(argv,index,defval)
%
	if (length(argv)>=index)
		val=argv(){index};
		if (ischar(defval))
			ret=val;
		else
			ret=str2num(argv(){index});
		endif
	else
		ret=defval;
	endif
endfunction