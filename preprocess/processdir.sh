#!/bin/sh
find $1 -iname "p_tpd_*_??.asc" -exec preprocess.sh \{\} \;

if [ -d $1/EXPORT ]; then
  preprocessIR.m $1/IR.irdat $1/EXPORT/
fi

ls -l $1/*.dat
