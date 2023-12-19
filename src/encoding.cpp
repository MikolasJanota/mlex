/*
 * File:  encoding.cpp
 * Author:  mikolas
 * Created on:  Tue Dec 13 12:20:37 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "encoding.h"
#include "auxiliary.h"
#include "binary_function.h" // for BinaryFunction
#include "minisat/core/SolverTypes.h"
#include <cassert> // for assert
#include <vector>
using SATSPC::lit_Undef;

inline bool is_true(SATSPC::Lit l, const SATSPC::vec<SATSPC::lbool> &values) {
    const auto v = SATSPC::var(l);
    if (v >= values.size())
        return false;
    const auto value = values[v];
    return !SATSPC::sign(l) ? (value == SATSPC::l_True)
                            : (value == SATSPC::l_False);
}

void Encoding::encode_bij() {
    const auto n = d_table.order();
    std::vector<SATSPC::Lit> ls;
#ifdef USE_MINISATSIMP
    for (size_t d = 0; d < n; d++)
        for (size_t r = 0; r < n; r++)
            d_sat.setFrozen(var(perm(d, r)), true);
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
}

void Encoding::encode_shot(const std::pair<size_t, size_t> &cell,
                           const std::vector<size_t> &vals,
                           SATSPC::Lit selector) {
    assert(selector != lit_Undef);
    selector = ~selector;
    const auto n = d_table.order();
    const auto &[row, col] = cell;
    auto &ls = _encoding_pos_ls;
    if (row == col) {
        for (size_t e = 0; e < n; e++) {
            if (!d_table.is_set(e, e))
                continue;
            ls.clear();
            ls.push(~perm(e, row));
            const auto old_val = d_table.get(e, e);
            bool taut = false;
            for (const auto &val : vals) {
                if (old_val == e && row != val)
                    continue;
                if (old_val == e && row == val) {
                    taut = true;
                    break;
                }
                ls.push(perm(old_val, val));
            }
            if (taut)
                continue;
            ls.push(selector);
            d_sat.addClause_(ls);
        }
    } else {
        for (size_t dom_row = 0; dom_row < n; dom_row++)
            for (size_t dom_col = 0; dom_col < n; dom_col++) {
                if (dom_col == dom_row)
                    continue;
                if (!d_table.is_set(dom_row, dom_col))
                    continue;
                ls.clear();
                ls.push(~perm(dom_row, row));
                ls.push(~perm(dom_col, col));
                const auto old_val = d_table.get(dom_row, dom_col);
                bool taut = false;
                for (const auto &val : vals) {
                    if ((old_val == dom_row && row != val) ||
                        (old_val == dom_col && col != val))
                        continue;
                    if ((old_val == dom_row && row == val) ||
                        (old_val == dom_col && col == val)) {
                        taut = true;
                        break;
                    }
                    ls.push(perm(old_val, val));
                }
                if (taut)
                    continue;
                ls.push(selector);
                d_sat.addClause_(ls);
            }
    }
}

void Encoding::encode_pos(const Assignment &assignment, SATSPC::Lit selector) {
    const auto &[row, col, val] = assignment;
    const auto n = d_table.order();
    const auto sel = selector != lit_Undef;
    if (sel)
        selector = ~selector;
    auto &ls = _encoding_pos_ls;
    if (row == col) {
        for (size_t e = 0; e < n; e++) {
            if (!d_table.is_set(e, e))
                continue;
            const auto old_val = d_table.get(e, e);
            if (old_val == e && row == val)
                continue; // tautology
            ls.clear();
            ls.push(~perm(e, row));
            const auto contradiciton = old_val == e && val != row;
            if (!contradiciton)
                ls.push(perm(old_val, val));
            if (selector != lit_Undef)
                ls.push(selector);
            d_sat.addClause_(ls);
        }
    } else {
        for (size_t dom_row = 0; dom_row < n; dom_row++)
            for (size_t dom_col = 0; dom_col < n; dom_col++) {
                if (dom_col == dom_row)
                    continue;
                if (!d_table.is_set(dom_row, dom_col))
                    continue;
                const auto old_val = d_table.get(dom_row, dom_col);
                if ((old_val == dom_row && val == row) ||
                    (old_val == dom_col && val == col))
                    continue; // tautology
                const auto contradiciton = (old_val == dom_row && val != row) ||
                                           (old_val == dom_col && val != col);
                ls.clear();
                int sz = 2 + (sel ? 1 : 0) + (contradiciton ? 0 : 1);
                ls.growTo(sz);
                if (selector != lit_Undef)
                    ls[--sz] = selector;
                if (!contradiciton)
                    ls[--sz] = perm(old_val, val);
                ls[--sz] = ~perm(dom_col, col);
                ls[--sz] = ~perm(dom_row, row);
                assert(sz == 0);
                d_sat.addClause_(ls);
            }
    }
}

void Encoding::print_solution(std::ostream &output) {
    const auto n = d_table.order();
    const auto &model = d_sat.model();
    std::vector<size_t> p(n, 0);
    for (size_t a = 0; a < n; a++)
        for (size_t b = 0; b < n; b++)
            if (is_true(perm(a, b), model)) {
                p[a] = b;
                break;
            }
    output << "[ " << std::endl;
    for (size_t row = 0; row < n; row++) {
        output << "[ ";
        for (size_t col = 0; col < n; col++) {
            if (col)
                output << " , ";
            output << p[d_table.get(row, col)] + 1;
        }
        output << "]" << std::endl;
        ;
    }
    output << "]" << std::endl;
}

