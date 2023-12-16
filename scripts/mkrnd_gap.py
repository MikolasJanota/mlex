#!/usr/bin/env python3
# File:  mkrnd.py
# Author:  mikolas
# Created on:  Mon Aug 14 11:38:40 AM UTC 2023
# Copyright (C) 2023, Mikolas Janota

import sys
import random

def prn(g, f):
    print('[', end='', file=f)
    ord=len(g)
    for rx,r in enumerate(g):
        print(f"{' ' if rx > 0 else ''}[", end='', file=f)
        for cx,c in enumerate(r):
            print(f"{',' if cx > 0 else ''}{c+1}", end='', file=f)
        print(']', end='', file=f)
        if rx+1<ord:
           print(',', file=f)
    print(']', file=f, end='')

def run_main():
    """ run the whole thing """
    if len(sys.argv) != 4:
        print("Expected ORDER COUNT RANDOM-SEED")
        return 1
    order = int(sys.argv[1])
    count = int(sys.argv[2])
    seed = int(sys.argv[3])
    random.seed(seed)
    rng = range(order)

    with open(f"rnd_{order}_{count}_{seed}.gap", 'w') as f:
        print('[', file=f)
        for c in range(count):
            prn([[random.randint(0,order-1) for col in rng] for row in rng],f)
            if c+1 < count:
              print(',', file=f)
        print('\n]', file=f)

if __name__ == "__main__":
    sys.exit(run_main())
