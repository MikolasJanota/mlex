#!/usr/bin/env python3
# File:  bruteforce.py
# Author:  mikolas
# Created on:  Sat Dec 16 19:25:47 CET 2023
# Copyright (C) 2023, Mikolas Janota
from itertools import permutations
from statistics import mean
import ast
import sys
import copy

def inv(pi):
    """ Calculate the inverse of a given permutation. """
    rv = [-1] * len(pi)
    for i, v in enumerate(pi):
        rv[v] = i
    return rv

def prn(g):
    """ Pretty-print a binary table. """
    print('[', end='')
    ord=len(g)
    for rx,r in enumerate(g):
        print(f"{' ' if rx > 0 else ''}[", end='')
        for cx,c in enumerate(r):
            print(f"{',' if cx > 0 else ''}{c+1}", end='')
        print(']', end='')
        if rx+1<ord:
           print(',')
    print(']', end='')

def perm(g, pi, pi1):
    """ Calculate the isomorphic copy of g given a permutation and its inverse. """
    ord = len(g)
    rng = range(ord)
    rv = [[-1] * ord for _ in rng]
    for r in rng:
        for c in rng:
            rv[r][c] = pi[g[pi1[r]][pi1[c]]]
    return rv

def cmp(orig, ming, pi, pi1):
    """ Compare the isomorphic copy of orig to ming. """
    rng = range(len(orig))
    for r in rng:
        for c in rng:
            val = ming[r][c]
            nval = pi[orig[pi1[r]][pi1[c]]]
            if nval < val:
                return -1
            if nval > val:
                return +1
    return 0

def prn_perm(pi):
    """ Pretty-print a permutation. """
    ord = len(pi)
    seen = set()
    rv = ''
    for i in range(ord):
        if i in seen:
            continue
        cy = []
        v = i
        while v not in seen:
            seen.add(v)
            cy.append(v+1)
            v = pi[v]
        rv += '(' + ' '.join(map(str, cy)) + ')'
    return rv

def solve(g):
    """ Find the smallest isomorphic copy of g. """
    rng = range(len(g))
    m = copy.deepcopy(g)
    best_pi = [i for i in rng]
    updates, last = 0,0
    for i, pi in enumerate(permutations(rng)):
        pi1 = inv(pi)
        if cmp(g, m, pi, pi1) < 0:
            last = i
            m = perm(g, pi, pi1)
            best_pi = copy.deepcopy(pi)
            updates += 1
    prn(m)
    return updates, last, best_pi

def run_main():
    updates_list = []
    last_list = []
    gs = ast.literal_eval(sys.stdin.read())
    print('[')
    for i,g in enumerate(gs):
        for rx,r in enumerate(g):
            for cx,c in enumerate(r):
               r[cx]-=1
        updates, last, pi = solve(g)
        if i+1 < len(gs):
           print(',', flush=True)
        else:
           print(flush=True)
        print(f"# updates: {updates}", flush=True)
        print(f"# pi {last}: {prn_perm(pi)}")
        updates_list.append(updates)
        last_list.append(last)
    print(']')
    print(f"# avg updates: {mean(updates_list)}")
    print(f"# avg permutation: {mean(last_list)}")

if __name__ == "__main__":
    sys.exit(run_main())
