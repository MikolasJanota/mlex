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

set -e
INPUT_FILE=$1
echo '% print a nauty file for each model'
./mlex -G -m $INPUT_FILE
echo '% run nauty on all the files'
ls *.dre | parallel -n1 ./rd.sh 
echo '% move outputs to a new directory'
mkdir outs
mv *.dre.out outs/
echo '% look for duplicate files'
cd outs/
fdupes -q -f . >out_dup
echo '% DONE'
DUPS=`grep -vc '^$' out_dup`
ALL=`ls *.dre.out| wc -l`
NONISO=$((ALL-DUPS))
echo "Getting $NONISO noniso out of $ALL"
