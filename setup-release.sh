#!/bin/bash
#
# File:  setup-release.sh
# Author:  mikolas
# Created on:  Sat Jul 31 20:17:28 WEST 2021
# Copyright (C) 2021, Mikolas Janota
#
set -e
set -u
if [ ! -d "cadical" ]; then
  echo setting up a SAT solver
  ./setup-cadical.sh
  echo SAT solver done
fi

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
echo "all done, binary in build"
