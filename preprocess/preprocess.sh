#!/bin/sh

PRESS=$1
TEMP=/tmp
TPDTools=~/octave/TPDTools

OUTFILE=`echo $1 | sed -e "s/.asc/.dat/g" | sed -e "s/p_tpd/tpd/g"`
DOSE=`echo $1 | sed -e "s/p_tpd_/q_dose_/g"`
DOSEP=`echo $1 | sed -e "s/p_tpd_/p_dose_/g"`
TPD=`echo $1 | sed -e "s/p_tpd_/q_tpd_/g"`

echo "writing to " $OUTFILE

cat $PRESS | head -n 1 > $TEMP/time.txt
cat $PRESS | sed -e "s/,/./g" | tail -n+2 > $TEMP/press.txt
cat $TPD | sed -e "s/,/./g"  > $TEMP/tpd.txt
cat $DOSE | sed -e "s/,/./g"  > $TEMP/dose.txt
cat $DOSEP | sed -e "s/,/./g" | tail -n+2 > $TEMP/dosep.txt

octave -q $TPDTools/preprocess/preprocesslv.m $OUTFILE $TEMP
