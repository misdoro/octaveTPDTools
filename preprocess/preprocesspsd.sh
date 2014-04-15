#!/bin/sh

PRESS=$1
TEMP=/tmp
TPDTools=/home/doronin/octave/TPDTools

OUTFILE=`echo $1 | sed -e "s/.asc/.dat/g" | sed -e "s/p_psd/idle/g"`
TPD=`echo $1 | sed -e "s/p_psd_/q_psd_/g"`

echo "writing to " $OUTFILE

cat $PRESS | sed -e "s/,/./g" | tail -n+2 > $TEMP/press.txt
cat $TPD | sed -e "s/,/./g"  > $TEMP/tpd.txt

octave -q $TPDTools/preprocess/preprocesspsd.m $OUTFILE $TEMP
