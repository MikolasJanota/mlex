/*
 * File:  encoding_explicit.cpp
 * Author:  mikolas
 * Created on:  Tue Jul 18 16:54:16 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "encoding_explicit.h"
#include "minisat/core/SolverTypes.h"
#include <cassert> // for assert
#include <cstddef>
#include <vector>

inline bool is_true(SATSPC::Lit l, const SATSPC::vec<SATSPC::lbool> &values) {
    const auto v = SATSPC::var(l);
    if (v >= values.size())
        return false;
    const auto value = values[v];
    return !SATSPC::sign(l) ? (value == SATSPC::l_True)
                            : (value == SATSPC::l_False);
}

SATSPC::Lit EncodingExplicit::encode_less(const BinaryFunction &table) {
    return encode_less_rec(table, 0, 0);
}

void EncodingExplicit::make_solution(BinaryFunction &table) {
    const auto n = d_table.order();
    const auto &model = d_sat.model();
    assert(table.order() == n);
    for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++) {
            for (size_t v = 0; v < n; v++)
                if (is_true(val(r, c, v), model)) {
                    table.set(r, c, v);
                    break;
                }
            assert(table.is_set(r, c));
        }
}

SATSPC::Lit EncodingExplicit::encode_less_rec(const BinaryFunction &table,
                                              size_t r, size_t c) {
    const auto n = d_table.order();
    if (r == n)
        return ~d_sat.true_lit();
    const size_t nc = (c + 1) % n;
    const size_t nr = nc == 0 ? r + 1 : r;
    const auto nl = encode_less_rec(table, nr, nc);
    const auto l = lesslit();
    const auto v = table.get(r, c);
    // disable all values above the current one
    for (size_t i = v + 1; i < n; i++)
        d_sat.addClause(~l, ~val(r, c, i));

    d_sat.addClause(~l, ~val(r, c, v), nl);
    return l;
}

void EncodingExplicit::encode_iso() {
    const auto n = d_table.order();
    SATSPC::vec<SATSPC::Lit> ls;
    for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++) {
            const auto v = d_table.get(r, c);
            for (size_t r1 = 0; r1 < n; r1++)
                for (size_t c1 = 0; c1 < n; c1++)
                    for (size_t v1 = 0; v1 < n; v1++) {
                        if ((r == c && r1 != c1) || (r == v && r1 != v1) ||
                            (c == v && c1 != v1))
                            continue;
                        if ((r1 == c1 && r != c) || (r1 == v1 && r != v) ||
                            (c1 == v1 && c != v))
                            continue;
                        ls.clear();
                        ls.push(~perm(r, r1));
                        ls.push(~perm(c, c1));
                        ls.push(~perm(v, v1));
                        ls.push(val(r1, c1, v1));
                        d_sat.addClause_(ls);
                    }
        }
}

void EncodingExplicit::encode_bij() {
    const auto n = d_table.order();
    std::vector<SATSPC::Lit> ls;
#ifdef USE_MINISATSIMP
    for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++)
            for (size_t v = 0; v < n; v++)
                d_sat.setFrozen(var(val(r, c, val)), true);
#endif

    for (size_t d = 0; d < n; d++) {
        ls.clear();
        for (size_t r = 0; r < n; r++)
            ls.push_back(perm(d, r));
        eq1(d_sat, ls);
    }

    for (size_t r = 0; r < n; r++) {
        ls.clear();
        for (size_t d = 0; d < n; d++)
            ls.push_back(perm(d, r));
        eq1(d_sat, ls);
    }

    for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++) {
            ls.clear();
            for (size_t v = 0; v < n; v++)
                ls.push_back(val(r, c, v));
            eq1(d_sat, ls);
        }
}

