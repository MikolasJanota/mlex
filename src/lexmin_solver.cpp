/*
 * File:  lexmin_solver.cpp
 * Author:  mikolas
 * Created on:  Wed Dec 14 13:29:54 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "lexmin_solver.h"
#include "auxiliary.h"
#include "minisat/core/SolverTypes.h"
#include "minisat_ext.h"
#include <cassert>
#include <math.h>
#include <memory>
#include <vector>
using SATSPC::Lit;
using SATSPC::mkLit;

#ifdef NDEBUG
#define TRACE(code)
#else
#define TRACE(code)                                                            \
    do {                                                                       \
        code                                                                   \
    } while (0)
#endif

void LexminSolver::solve() {
    if (d_options.incremental) {
        d_sat      = std::make_unique<SATSPC::MiniSatExt>();
        d_encoding = std::make_unique<Encoding>(d_options, *d_sat, d_table);
        d_encoding->encode_bij();
    }

    const auto n = d_table.order();

    d_assignments.clear();
    d_assignments.reserve(n * n);

    for (size_t row = 0; row < n; row++) {
        for (size_t col = 0; col < n; col++) {
            d_assignments.push_back({row, col, 0});
            TRACE(d_output.comment(3) << row << " " << col << " :";);
            auto &current_assignment = d_assignments.back();
            auto &cur_val            = std::get<2>(current_assignment);
            while (!test_sat()) {
                cur_val++;
                assert(cur_val < n);
            }
            TRACE(d_output.ccomment(3) << std::endl;);
        }
    }
    make_solution();
}

bool LexminSolver::test_sat() {
    const auto start_time = read_cpu_time();
    const auto rv = d_options.incremental ? test_sat_inc() : test_sat_noinc();
    d_statistics.satTime->inc(read_cpu_time() - start_time);
    d_statistics.satCalls->inc();
    TRACE(d_output.ccomment(3)
              << " " << std::get<2>(d_assignments.back()) << ":"
              << SHOW_TIME(read_cpu_time() - start_time););
    return rv;
}

bool LexminSolver::test_sat_inc() {
    Minisat::vec<Minisat::Lit> assumps(1, mkLit(d_sat->fresh()));
    const auto &               selector = assumps[0];
    d_encoding->encode_pos(d_assignments.back(), selector);
    const auto res = d_sat->solve(assumps);
    d_sat->addClause(res ? selector : ~selector);
    return res;
}

bool LexminSolver::test_sat_noinc() {
    d_sat      = std::make_unique<SATSPC::MiniSatExt>();
    d_encoding = std::make_unique<Encoding>(d_options, *d_sat, d_table);
    d_encoding->encode(d_assignments);
    const auto res = d_sat->solve();
    d_encoding.reset();
    d_sat.reset();
    return res;
}

void LexminSolver::make_solution() {
    const auto n = d_table.order();
    d_solution   = std::make_unique<BinaryFunction>(n);
    d_solution->set_name(d_table.get_name());
    d_solution->set_additional_info(d_table.get_additional_info());
    for (const auto &[row, col, val] : d_assignments)
        d_solution->set(row, col, val);
}

void LexminSolver::print_solution(std::ostream &output) {
    if (d_options.mace_format)
        d_solution->print_mace(output);
    else
        d_solution->print_gap(output);
}

