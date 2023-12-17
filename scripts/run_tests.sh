#!/bin/bash
#
# File:  run_tests.sh
# Author:  mikolas
# Created on:  Sat Dec 16 09:11:29 PM UTC 2023
# Copyright (C) 2023, Mikolas Janota
#
#
#
SOLVER=../build/mlex

if ! ( which `cut <<<${SOLVER} -f1 -d\ ` >/dev/null 2>&1 ) ; then
    echo "Looks like $SOLVER can't be run."
    exit 300
fi

TMPF=/tmp/${RANDOM}.gap
TMPF1=/tmp/${RANDOM}.gap
ls ../tests/*.gap.gz | grep -v '_min\.gap.gz' |\
  while read F;  do
    echo $F
    $SOLVER $F | ./reformat.py >$TMPF
    zcat ${F%.gap.gz}_min.gap.gz >$TMPF1
    if diff $TMPF1 $TMPF; then
      echo okay
    else
      echo KO
    fi
  done
rm -f $TMPF $TMPF1
