#!/usr/bin/env python3
""" Small utility to reformat tables. """

# File:  reformat.py
# Author:  mikolas
# Created on:  Sat Dec 16 21:20:49 CET 2023
# Copyright (C) 2023, Mikolas Janota
import sys
import ast


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


def run_main():
    """Run the whole program."""
    gs = ast.literal_eval(sys.stdin.read())
    print("[")
    for i, g in enumerate(gs):
        prn(g)
        if i + 1 < len(gs):
            print(",")
    print("\n]")
    return 0


if __name__ == "__main__":
    sys.exit(run_main())
