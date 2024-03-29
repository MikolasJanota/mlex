#!/usr/bin/env python3
""" Brute force solver for the lexmin problem. """

# File:  bruteforce.py
# Author:  mikolas
# Created on:  Sat Dec 16 19:25:47 CET 2023
# Copyright (C) 2023, Mikolas Janota
from itertools import permutations
from itertools import product
from statistics import mean
import argparse
import ast
import sys
import copy


def inv(pi):
    """Calculate the inverse of a given permutation."""
    rv = [-1] * len(pi)
    for i, v in enumerate(pi):
        rv[v] = i
    return rv


def prn(g):
    """Pretty-print a binary table."""
    print("[", end="")
    order = len(g)
    for rx, r in enumerate(g):
        print(f"{' ' if rx > 0 else ''}[", end="")
        for cx, c in enumerate(r):
            print(f"{',' if cx > 0 else ''}{c+1}", end="")
        print("]", end="")
        if rx + 1 < order:
            print(",")
    print("]", end="")


def perm(g, pi, pi1):
    """Calculate the isomorphic copy of g given a permutation and its inverse."""
    order = len(g)
    rng = range(order)
    rv = [[-1] * order for _ in rng]
    for r in rng:
        for c in rng:
            rv[r][c] = pi[g[pi1[r]][pi1[c]]]
    return rv


def cmp(orig, ming, pi, pi1):
    """Compare the isomorphic copy of orig to ming."""
    rng = range(len(orig))
    for r, c in product(rng, rng):
        val = ming[r][c]
        nval = pi[orig[pi1[r]][pi1[c]]]
        if nval < val:
            return -1
        if nval > val:
            return +1
    return 0


def prn_perm(pi):
    """Pretty-print a permutation."""
    order = len(pi)
    seen = set()
    rv = ""
    for i in range(order):
        if i in seen:
            continue
        cy = []
        v = i
        while v not in seen:
            seen.add(v)
            cy.append(v + 1)
            v = pi[v]
        rv += "(" + " ".join(map(str, cy)) + ")"
    return rv


def count_equal(g, ming):
    """Count how many permutations turn g into ming."""
    assert len(g) == len(ming)
    rng = range(len(g))
    rv = 0
    for pi in permutations(rng):
        pi1 = inv(pi)
        if cmp(g, ming, pi, pi1) == 0:
            rv += 1
    return rv


def solve(args, g):
    """Find the smallest isomorphic copy of g."""
    order = len(g)
    rng = range(order)
    m = copy.deepcopy(g)
    best_pi = list(rng)
    updates, last = 0, 0
    for i, pi in enumerate(permutations(rng)):
        pi1 = inv(pi)
        if cmp(g, m, pi, pi1) < 0:
            last = i
            m = perm(g, pi, pi1)
            best_pi = copy.deepcopy(pi)
            updates += 1
    prn(m)
    solutions = count_equal(g, m) if args.count else -1
    return updates, last, best_pi, solutions


def solve_all(args, gs):
    """Solve all the tables gs."""
    updates_list = []
    last_list = []
    print("[")
    for i, g in enumerate(gs):
        order = len(g)
        print(f"# order: {order}")
        for r in g:
            for cx, _ in enumerate(r):
                r[cx] -= 1
                assert 0 <= r[cx] < order
        updates, last, pi, solutions = solve(args, g)
        if i + 1 < len(gs):
            print(",", flush=True)
        else:
            print(flush=True)
        print(f"# updates: {updates}", flush=True)
        print(f"# pi {last}: {prn_perm(pi)}")
        if args.count:
            print(f"# solutions: {solutions}")
        updates_list.append(updates)
        last_list.append(last)
    print("]")
    print(f"# avg updates: {mean(updates_list)}")
    print(f"# avg permutation: {mean(last_list)}")


def run_main():
    """Run the whole program."""
    arg_parser = argparse.ArgumentParser(
        description="Bruteforce solution for the minlex problem."
    )
    arg_parser.add_argument("-c", "--count", default=False, action="store_true")
    arg_parser.add_argument("filename", default="-", nargs="?")
    args = arg_parser.parse_args()

    if args.filename == "-":
        data = sys.stdin.read()
    else:
        with open(args.filename, "r", encoding="utf-8") as inf:
            data = inf.read()
    gs = ast.literal_eval(data)
    solve_all(args, gs)
    return 0


if __name__ == "__main__":
    sys.exit(run_main())
