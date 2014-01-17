#!/bin/sh
find $1 -iname "p_tpd_*_??.asc" -exec preprocess.sh \{\} \;
ls -l $1/*.dat
