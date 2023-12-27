#!/usr/bin/env python3
"""Generate random input tables."""

# File:  mkrnd.py
# Author:  mikolas
# Created on:  Mon Aug 14 11:38:40 AM UTC 2023
# Copyright (C) 2023, Mikolas Janota

import sys
import random


def run_main():
    """run the whole thing"""
    if len(sys.argv) != 4:
        print("Expected ORDER COUNT RANDOM-SEED")
        return 1
    order = int(sys.argv[1])
    count = int(sys.argv[2])
    seed = int(sys.argv[3])
    random.seed(seed)

    for c in range(count):
        with open(f"rnd_{seed}_{order}_{c}.out", encoding="utf-8", mode="w") as f:
            f.write(f"interpretation({order}, [number={count}, seed={seed}], [\n")
            f.write(" function(*(_,_), [\n")
            for i in range(order):
                f.write(
                    3 * " "
                    + ",".join(
                        map(str, [random.randint(0, order - 1) for j in range(order)])
                    )
                )
                if i + 1 < order:
                    f.write(",\n")
            f.write("])]).\n")

    return 0


if __name__ == "__main__":
    sys.exit(run_main())
