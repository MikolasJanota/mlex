#!/usr/bin/env python3
# File:  bruteforce.py
# Author:  mikolas
# Created on:  Sat Dec 16 19:25:47 CET 2023
# Copyright (C) 2023, Mikolas Janota
from itertools import permutations
import ast
import sys
import copy
update_list = []

def inv(pi):
    rv = [-1] * len(pi)
    for i, v in enumerate(pi):
        rv[v] = i
    return rv

def prn(g):
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
    ord = len(g)
    rng = range(ord)
    rv = [[-1] * ord for _ in rng]
    for r in rng:
        for c in rng:
            rv[r][c] = pi[g[pi1[r]][pi1[c]]]
    return rv

def cmp(orig, ming, pi, pi1):
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
    rng = range(len(g))
    m = copy.deepcopy(g)
    best_pi = [i for i in rng]
    updates = 0
    for pi in permutations(rng):
        pi1 = inv(pi)
        if cmp(g, m, pi, pi1) < 0:
            m = perm(g, pi, pi1)
            best_pi = copy.deepcopy(pi)
            updates += 1
    prn(m)
    update_list.append(updates)
    return best_pi

def run_main():
    gs = ast.literal_eval(sys.stdin.read())
    print('[')
    for i,g in enumerate(gs):
        for rx,r in enumerate(g):
            for cx,c in enumerate(r):
               r[cx]-=1
        pi = solve(g)
        if i+1 < len(gs):
           print(',', flush=True)
        else:
           print(flush=True)
        print(f"# updates: {update_list[-1]}", file=sys.stderr, flush=True)
        print(f"# pi: {prn_perm(pi)}")
    print(']')
    print(f"avg updates: {sum(update_list)/len(update_list)}", file=sys.stderr)

if __name__ == "__main__":
    sys.exit(run_main())
