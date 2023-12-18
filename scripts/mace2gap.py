#!/usr/bin/env python3
# File:  mace2gap.py
# Author:  mikolas
# Created on:  Sun Dec 17 22:37:48 CET 2023
# Copyright (C) 2023, Mikolas Janota
import sys
import ast

def prn(g):
    print('[', end='')
    ord=len(g)
    for rx,r in enumerate(g):
        print(f"{' ' if rx > 0 else ''}[", end='')
        for cx,c in enumerate(r):
            print(f"{',' if cx > 0 else ''}{c}", end='')
        print(']', end='')
        if rx+1<ord:
           print(',')
    print(']', end='')

def read():
    gs = []
    for l in sys.stdin:
        l = l.strip()
        if not l:
            continue
        if l.startswith('interpretation'):
            continue
        if l.startswith('function'):
            g = []
            ord = None
            continue
        last = l.endswith(' ])]).')
        if last:
            l = l.split()[0].strip()
        else:
            assert l.endswith(',')
            l = l[:-1]

        row = [x+1 for x in map(int,l.split(','))]
        if ord is None:
            ord = len(row)
        else:
            assert len(row) == ord
        assert all(map(lambda x: 1<=x<=ord, row))
        g.append(row)
        if last:
            assert len(g) == ord
            gs.append(g)
            g = None
    return gs


def run_main():
    gs = read()
    print('[')
    for i,g in enumerate(gs):
        prn(g)
        if i+1 < len(gs):
           print(',', flush=True)
    print('\n]', flush=True)

if __name__ == "__main__":
    sys.exit(run_main())
