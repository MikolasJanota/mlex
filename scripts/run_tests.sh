#!/bin/bash
#
# File:  run_tests.sh
# Author:  mikolas
# Created on:  Sat Dec 16 09:11:29 PM UTC 2023
# Copyright (C) 2023, Mikolas Janota
#
#

TMPF=/tmp/${RANDOM}.gap
ls ../tests/*.gap | grep -v '_min\.gap' |\
  while read F;  do
    echo $F
    ../build/mlex $F | ./reformat.py >$TMPF
    if diff ${F%.gap}_min.gap $TMPF; then
      echo okay
    else
      echo KO
    fi
  done
rm -f $TMPF
