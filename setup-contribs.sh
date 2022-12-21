#!/bin/bash
#
# File:  setup-contribs.sh
# Author:  mikolas
# Created on:  Wed Dec 21 17:36:13 CET 2022
# Copyright (C) 2022, Mikolas Janota
#

set -e
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
