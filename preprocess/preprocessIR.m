#!/usr/bin/octave -q
if (nargin>=1);
	irpath=argv(){1};
  outfile="IR.irdat";
else
  printf("Usage: preprocessIR /path/to/irdata/\n");
  exit();
endif;

IRDAT=loadIRData(irpath);
printf("Loaded %d IR data files\n",IRDAT.count);

printf("Saving IR data to %s\n",outfile);
save("-binary",outfile,"IRDAT");