#!/usr/bin/octave -q
if (nargin>=2);
	irpath=argv(){2};
  outfile=argv(){1};
else
  printf("Usage: preprocessIR out.irdat /path/to/irdata/\n");
  exit();
endif;

IRDAT=loadIRData(irpath);
printf("Loaded %d IR data files\n",IRDAT.count);

printf("Saving IR data to %s\n",outfile);
save("-binary",outfile,"IRDAT");