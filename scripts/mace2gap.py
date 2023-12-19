#!/usr/bin/env python3
""" Small utility to convert MACE format to GAP format. """

# File:  mace2gap.py
# Author:  mikolas
# Created on:  Sun Dec 17 22:37:48 CET 2023
# Copyright (C) 2023, Mikolas Janota
import sys


def prn(g):
    """Pretty-print a binary table."""
    print("[", end="")
    order = len(g)
    for rx, r in enumerate(g):
        print(f"{' ' if rx > 0 else ''}[", end="")
        for cx, c in enumerate(r):
            print(f"{',' if cx > 0 else ''}{c}", end="")
        print("]", end="")
        if rx + 1 < order:
            print(",")
    print("]", end="")


def read():
    """Read tables on stdin."""
    gs = []
    for ln in sys.stdin:
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith("interpretation"):
            continue
        if ln.startswith("function"):
            g = []
            order = None
            continue
        last = ln.endswith(" ])]).")
        if last:
            ln = ln.split()[0].strip()
        else:
            assert ln.endswith(",")
            ln = ln[:-1]

        row = [x + 1 for x in map(int, ln.split(","))]
        if order is None:
            order = len(row)
        assert len(row) == order
        assert all(map(lambda x: 1 <= x <= order, row))
        g.append(row)
        if last:
            assert len(g) == order
            gs.append(g)
            g = None
    return gs


def run_main():
    """Run the whole program."""
    gs = read()
    print("[")
    for i, g in enumerate(gs):
        prn(g)
        if i + 1 < len(gs):
            print(",", flush=True)
    print("\n]", flush=True)


if __name__ == "__main__":
    sys.exit(run_main())
