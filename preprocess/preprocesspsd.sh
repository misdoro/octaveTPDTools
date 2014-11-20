#!/bin/sh
ARGC=$#

if [ $ARGC -ge 2 ]; then
	PRESS=$1
	PREFIX=$2
	TEMP=/tmp
	TPDTools=~/octave/TPDTools

	OUTFILE=`echo $1 | sed -e "s/.asc/.dat/g" | sed -e "s/p_$PREFIX/$PREFIX/g"`
	TPD=`echo $1 | sed -e "s/p_${PREFIX}_/q_${PREFIX}_/g"`
	
	echo "writing to " $OUTFILE

	cat $PRESS | sed -e "s/,/./g" | tail -n+2 > $TEMP/press.txt
	cat $TPD | sed -e "s/,/./g"  > $TEMP/tpd.txt
	
	octave -q $TPDTools/preprocess/preprocesspsd.m $OUTFILE $TEMP

else
	echo "Usage: preprocesspsd.sh p_prefix_date_index.asc prefix"
fi
