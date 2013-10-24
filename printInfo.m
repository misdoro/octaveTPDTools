function printInfo(input);
	for idx=1:rows(input.sorted);
		filename=input.filenames{input.sorted(idx,1)};
		printf("\n-----------------\n%s\n",filename);
		load(filename);
		if (isfield(tpd,"iN"))
			printf("Temperature range: %3.1f to %3.1f K\n",min(tpd.T),max(tpd.T));
			
			rkm=tpd.rate*60;
			if (abs(rkm-round(rkm))>0.05)
				printf("Ramp rate: %4.2f K/min\n",tpd.rate*60);
			else
				printf("Ramp rate: %d K/min\n",round(tpd.rate*60));
			endif
			
			printf("Mid: %d\n",tpd.mids(find(tpd.mids>0)));
		else
			printf("File is in old or unknown format\n")
		endif
	end
endfunction