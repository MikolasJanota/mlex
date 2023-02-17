#!/bin/bash
#
# File:  cdre.sh
# Author:  mikolas
# Created on:  Thu Feb 16 21:26:21 CET 2023
# Copyright (C) 2023, Mikolas Janota
#
#
#  The script rd.sh needs to be in the current directory and mlex dreadnaut as well.
#  fdupes need to be installed which might not always be the default
#  (or just symlinked)
#
#
#  The input file should also be in the current directory otherwise outputs
#  are created in the other directory, which descriptors and expect
set -u

MLEX_OPTS="-m"
INPUT_FILE=$1
echo '% print a nauty file for each model'
./mlex ${MLEX_OPTS} -G $INPUT_FILE
echo '% run nauty on all the files'
ls *.dre | parallel -n1 ./rd.sh 
OUTPUTS=outs${RANDOM}
echo "% move outputs to a new directory ${OUTPUTS}"
mkdir ${OUTPUTS}
mv *.dre.out ${OUTPUTS}
echo '% look for duplicate files'
COUNT_ALL=`ls ${OUTPUTS}/*.dre.out | wc -l`
let COUNT_DUPS=`grep -vc '^$' out_dup`
fdupes -q -f ${OUTPUTS} >out_dup
COUNT_NONISO=$((COUNT_ALL-COUNT_DUPS))
echo "% done"
echo "Getting ${COUNT_NONISO} noniso out of ${COUNT_ALL}."
