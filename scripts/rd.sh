#!/bin/bash
#
# File:  rd.sh
# Author:  mikolas
# Created on:  Sun Feb 12 10:26:34 WET 2023
# Copyright (C) 2023, Mikolas Janota
#
f=$1
(time ./dreadnaut <${f} 2>${f}.err | grep -e '[[:digit:]]\+[[:space:]]*:.*;' >${f}.out ) 2>${f}.time
