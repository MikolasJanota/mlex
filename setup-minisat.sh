#!/bin/bash
#
# File:  setup-minisat.sh
# Author:  mikolas
# Created on:  Wed Dec 21 17:13:24 CET 2022
# Copyright (C) 2022, Mikolas Janota
#

set -e
mkdir -p contrib
cd contrib
mkdir minisat 
cd minisat 
git clone git@github.com:agurfinkel/minisat.git
cd minisat 
make config prefix=..
make -j4 install
