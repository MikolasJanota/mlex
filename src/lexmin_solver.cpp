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

void LexminSolver::solve() {
    if (d_options.incremental) {
        d_sat      = std::make_unique<SATSPC::MiniSatExt>();
        d_encoding = std::make_unique<Encoding>(d_options, *d_sat, d_table);
        d_encoding->encode_bij();
    }

    const auto n = d_table.order();
    for (size_t row = 0; row < n; row++) {
        for (size_t col = 0; col < n; col++) {
            if (d_options.verbose > 1)
                comment() << row << " " << col << " :";
            bool found = false;
            for (size_t val = 0; !found && val < n; val++) {
                found = test_sat({row, col, val});
            }
            assert(found);
            if (d_options.verbose > 1)
                std::cout << std::endl;
        }
    }
    make_solution();
}

bool LexminSolver::test_sat(const Encoding::Assignment &assignment) {
    const auto start_time = read_cpu_time();
    const auto rv         = d_options.incremental ? test_sat_inc(assignment)
                                                  : test_sat_noinc(assignment);
    d_total_sat_time += read_cpu_time() - start_time;
    if (d_options.verbose > 1)
        std::cout << " " << std::get<2>(assignment) << ":"
                  << SHOW_TIME(read_cpu_time() - start_time);
    return rv;
}

bool LexminSolver::test_sat_inc(const Encoding::Assignment &assignment) {
    d_assignments.push_back(assignment);
    Minisat::vec<Minisat::Lit> assumps(1, mkLit(d_sat->fresh()));
    const auto &               selector = assumps[0];
    d_encoding->encode_pos(assignment, selector);
    const auto res = d_sat->solve(assumps);
    if (!res)
        d_assignments.pop_back();
    d_sat->addClause(res ? selector : ~selector);
    return res;
}

bool LexminSolver::test_sat_noinc(const Encoding::Assignment &assignment) {
    d_assignments.push_back(assignment);
    d_sat      = std::make_unique<SATSPC::MiniSatExt>();
    d_encoding = std::make_unique<Encoding>(d_options, *d_sat, d_table);
    if (d_options.verbose > 4) {
        comment() << " A:";
        for (const auto &[a, b, c] : d_assignments)
            std::cout << " (" << a << " " << b << " " << c << ")";
        std::cout << std::endl;
    }
    d_encoding->encode(d_assignments);
    const auto res = d_sat->solve();
    if (!res)
        d_assignments.pop_back();
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
    comment() << "total SAT time " << SHOW_TIME(d_total_sat_time) << "\n";
    const auto                       n = d_table.order();
    std::vector<std::vector<size_t>> t(n);
    for (size_t row = 0; row < n; row++)
        t[row].resize(n, static_cast<size_t>(-1));

    for (const auto &[row, col, val] : d_assignments) {
        assert(t[row][col] == static_cast<size_t>(-1));
        t[row][col] = val;
    }
    if (d_options.mace_format)
        d_solution->print_mace(output);
    else
        d_solution->print_gap(output);
}

