#!/bin/bash
#
# File:  setup-release.sh
# Author:  mikolas
# Created on:  Sat Jul 31 20:17:28 WEST 2021
# Copyright (C) 2021, Mikolas Janota
#

./setup-contribs.sh

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j4
echo "all done, binary in build"
