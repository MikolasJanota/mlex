/*
 * File:  lexmin_solver_explicit.cpp
 * Author:  mikolas
 * Created on:  Tue Jul 18 18:02:15 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "lexmin_solver_explicit.h"
#include "binary_function.h"
#include "encoding_explicit.h"
#include <cassert>
#include <cstddef>
#include <vector>
using SATSPC::Lit;
#if !defined(NDEBUG) || defined(SOLVER_TRACING)
#define TRACE(code)                                                            \
    do {                                                                       \
        code                                                                   \
    } while (0);                                                               \
    do {                                                                       \
        std::cout.flush();                                                     \
    } while (0)
#else
#define TRACE(code)
#endif

void LexminSolverExplicit::solve() {
    assert(!d_options.diagonal);
    const auto n = d_table.order();
    d_sat = std::make_unique<SATSPC::MiniSatExt>();
    d_encoding = std::make_unique<EncodingExplicit>(d_output, *d_sat, d_table);
    const auto start_time = read_cpu_time();
    d_encoding->encode_bij();
    d_encoding->encode_iso();
    d_statistics.encodingTime->inc(read_cpu_time() - start_time);
    if (d_options.opt1stRow && !d_options.diagonal)
        opt1stRow();

    if (d_solution)
        delete d_solution;
    d_solution = new BinaryFunction(n);
    d_solution->set(d_table);
    Lit last = d_encoding->encode_less(d_table);
    Minisat::vec<Minisat::Lit> assumps(1, last);
    while (run_sat(assumps)) {
        d_sat->addClause(~last);
        d_encoding->make_solution(*d_solution);
        TRACE(d_solution->print(comment(3) << "sol:\n", "#  ") << std::endl;);
        last = d_encoding->encode_less(*d_solution);
        assumps[0] = last;
    }
    if (d_options.verbose > 1) {
        std::vector<size_t> perm;
        make_permutation(perm);
        show_permutation(comment(2) << "perm:", perm) << std::endl;
    }
}

bool LexminSolverExplicit::run_sat(Minisat::vec<Minisat::Lit> &assumps) {
    const auto start_time = read_cpu_time();
    const auto rv = d_sat->solve(assumps);
    d_statistics.satTime->inc(read_cpu_time() - start_time);
    d_statistics.satCalls->inc();
    comment(3) << "sc:" << d_statistics.satCalls->get() << ":" << rv
               << SHOW_TIME(read_cpu_time() - start_time) << '\n';
    return rv;
}

CompFunction LexminSolverExplicit::make_solution_comp() {
    const auto n = d_table.order();
    CompFunctionBuilder b(n, 2);
    for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++)
            b.push(d_solution->get(r, c));
    return b.make();
}

void LexminSolverExplicit::make_permutation(std::vector<size_t> &perm) {
    const auto n = d_table.order();
    perm.clear();
    perm.resize(n, -1);
    for (auto dom = n; dom--;) {
        bool found = false;
        for (auto rng = n; !found && rng--;)
            if (d_sat->eval_lit(d_encoding->perm(dom, rng)) == SATSPC::l_True) {
                perm[dom] = rng;
                found = true;
            }
        assert(found);
    }
}

void LexminSolverExplicit::opt1stRow() {
    const auto n = d_table.order();
    // look for idempotents
    std::vector<size_t> idems;
    for (auto i = n; i--;)
        if (d_table.get(i, i) == i)
            idems.push_back(i);

    if (idems.empty())
        return; // TODO

    TRACE(print_set(comment(3) << "idems ", idems) << std::endl;);

    // count repetitions f(r,y)=r
    std::vector<size_t> repeats(n, 0);
    size_t max_repeats = 0;
    for (const auto row : idems) {
        for (auto col = n; col--;)
            if (d_table.get(row, col) == row)
                repeats[row]++;
        if (repeats[row] > max_repeats)
            max_repeats = repeats[row];
    }
    comment(3) << "max_repeats " << max_repeats << std::endl;

    assert(max_repeats > 0);

    // place zeros at the beginning of the first row
    for (size_t i = 0; i < max_repeats; i++) {
        d_sat->addClause(d_encoding->val(0, i, 0));
    }

    size_t some_first_row = -1; // some row which can be 1st
    size_t count_can_be_first = 0;
    for (auto row = n; row--;) {
        const bool can_be_first = repeats[row] == max_repeats;
        if (can_be_first) {
            count_can_be_first++;
            some_first_row = row;
        } else {
            d_sat->addClause(~d_encoding->perm(row, 0));
        }
    }

    assert(count_can_be_first > 0);
    comment(3) << "count_can_be_first " << count_can_be_first << std::endl;

    // handle the case of unique first row
    if (count_can_be_first == 1) {
        d_statistics.unique1stRow->inc();
        comment(2) << some_first_row << " fixed to 0 (opt1stRow)" << std::endl;
    }
}

