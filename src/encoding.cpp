/*
 * File:  encoding.cpp
 * Author:  mikolas
 * Created on:  Tue Dec 13 12:20:37 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "encoding.h"
#include "auxiliary.h"
#include "minisat/core/SolverTypes.h"
#include <vector>

inline bool is_true(SATSPC::Lit l, const SATSPC::vec<SATSPC::lbool> &values) {
    const auto v = SATSPC::var(l);
    if (v >= values.size())
        return false;
    const auto value = values[v];
    return !SATSPC::sign(l) ? (value == SATSPC::l_True)
                            : (value == SATSPC::l_False);
}

void Encoding::encode_bij() {
    const auto               n = d_table.order();
    std::vector<SATSPC::Lit> ls;
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
void Encoding::encode(const std::vector<Assignment> &assignments) {
    encode_bij();
    for (const auto &assignment : assignments)
        encode_pos(assignment, SATSPC::lit_Undef);
}

using SATSPC::lit_Undef;

void Encoding::encode_pos(const Assignment &assignment, SATSPC::Lit selector) {
    const auto &[row, col, val] = assignment;
    const auto n                = d_table.order();
    if (selector != lit_Undef)
        selector = ~selector;
    SATSPC::vec<SATSPC::Lit> ls;
    if (row == col) {
        for (size_t e = 0; e < n; e++) {
            const auto old_val = d_table.get(e, e);
            ls.clear();
            ls.push(~perm(e, row));
            ls.push(perm(old_val, val));
            if (selector != lit_Undef)
                ls.push(selector);
            d_sat.addClause_(ls);
            /* d_sat.addClause(~perm(e, row), perm(old_val, val)); */
        }
    } else {
        for (size_t dom_row = 0; dom_row < n; dom_row++)
            for (size_t dom_col = 0; dom_col < n; dom_col++) {
                if (dom_col == dom_row)
                    continue;
                const auto old_val = d_table.get(dom_row, dom_col);
                ls.clear();
                ls.push(~perm(dom_row, row));
                ls.push(~perm(dom_col, col));
                ls.push(perm(old_val, val));
                if (selector != lit_Undef)
                    ls.push(selector);
                d_sat.addClause_(ls);
                /* d_sat.addClause(~perm(dom_row, row), ~perm(dom_col, col), */
                /*                 perm(old_val, val)); */
            }
    }
}

void Encoding::print_solution(std::ostream &output) {
    const auto          n     = d_table.order();
    const auto &        model = d_sat.model();
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
