#!/usr/bin/env python3
# File:  rmdups.py
# Author:  mikolas
# Created on:  Mon Mar 27 12:56:31 CEST 2023
# Copyright (C) 2023, Mikolas Janota
import sys
import os
ll=''
for l in sys.stdin:
    l=l.strip()
    if l and ll:
        print("removing dup", l)
        os.remove(l)
    ll=l
