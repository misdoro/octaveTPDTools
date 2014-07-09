function printInfo(input,param);
	for idx=1:rows(input.sorted);
		filename=input.filenames{input.sorted(idx,1)};
		printf("\n-----------------\n%s\n",filename);
		load(filename);
    if (isfield(tpd,"version")&& tpd.version>=20140709)
      printf("TPD start at: %s",asctime(tpd.time));
    endif
		if (isfield(tpd,"iN"))
			if (isfield(tpd,"version")&& tpd.version>=20140120)
				for midx=find(dose.mids>0)
					printf("Mass %d dose integral %3.1e\n",dose.mids(midx),dose.integral(midx))
				end
				printf("Dose temperature %3.1f K\n", dose.T)
			else
				printf("Dose %3.1e A*s, at %3.1f K\n",dose.integral,dose.T);
			endif;
      
      printf("Temperature range: %3.1f to %3.1f K\n",min(tpd.T),max(tpd.T));
			rkm=tpd.rate*60;
			if (abs(rkm-round(rkm))>0.05)
				printf("Ramp rate: %4.2f K/min\n",tpd.rate*60);
			else
				printf("Ramp rate: %d K/min\n",round(tpd.rate*60));
			endif
			if(isfield(tpd,"version")&& tpd.version>=20140612)
        for midx=find(tpd.mids>0)
          printf("Mass %d tpd integral %3.1e (%3.1f ML)\n",
            tpd.mids(midx),tpd.intg(midx),tpd.intg(midx)/(param.monolayer));
        endfor
      else
  			printf("Estimated coverage %3.1f ML\n",tpd.integral/param.monolayer)
        printf("Mid: %d\n",tpd.mids(find(tpd.mids>0)));
			endif
			
			if (isfield(tpd,"version"))
				printf("File format version %d\n",tpd.version);
			endif
		else
			printf("File is in old or unknown format\n")
		endif
	end
endfunction