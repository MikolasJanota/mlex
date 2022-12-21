#!/bin/bash
#
# File:  setup-release.sh
# Author:  mikolas
# Created on:  Sat Jul 31 20:17:28 WEST 2021
# Copyright (C) 2021, Mikolas Janota
#
set -e
set -u
mkdir -p contrib

if [ ! -d "contrib/cadical" ]; then
  echo setting up a cadical SAT solver
  ./setup-cadical.sh
  echo cadical SAT solver done
fi

if [ ! -d "contrib/minisat" ]; then
  echo setting up a minisat SAT solver
  ./setup-minisat.sh
  echo minisat SAT solver done
fi

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
echo "all done, binary in build"
