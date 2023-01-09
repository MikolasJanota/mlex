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
    const auto n = d_table.order();
    d_fixed.resize(n, n);
    d_used.resize(n, false);
    if (d_options.incremental) {
        d_sat = std::make_unique<SATSPC::MiniSatExt>();
        make_encoding();
        d_encoding->encode_bij();
    }
    if (d_options.budgeting) {
        calculate_budgets();
    }

    d_assignments.clear();
    d_assignments.reserve(n * n);

    std::vector<size_t> current_row_budget;
    for (size_t row = 0; row < n; row++) {
        if (d_options.budgeting) {
            assert(d_rowBudget.size() == n);
            if (row == 1 && is_fixed(0)) {
                // mark first row as used
                d_used[0] = true;
                calculate_budgets();
            }
            current_row_budget = d_rowBudget;
        }

        for (size_t col = 0; col < n; col++) {
            d_assignments.push_back({row, col, 0});
            TRACE(d_output.comment(3) << "(" << row << " " << col << ") :";);
            auto &current_assignment = d_assignments.back();
            auto &cur_val = std::get<2>(current_assignment);
            while (true) {
                assert(cur_val < n);
                const bool found =
                    d_options.budgeting
                        ? current_row_budget[cur_val] && test_sat()
                        : test_sat();
                if (d_options.budgeting && found)
                    current_row_budget[cur_val]--;

                if (found)
                    break;
                cur_val++;
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

void LexminSolver::make_encoding() {
    d_encoding = std::make_unique<Encoding>(d_output, *d_sat, d_table);
    if (d_options.opt1stRow)
        opt1stRow();
}

void LexminSolver::opt1stRow() {
    // look for idempotents
    std::vector<size_t> idems;
    for (auto i = d_table.order(); i--;)
        if (d_table.get(i, i) == i)
            idems.push_back(i);

    // only idempotents can be at 1st row, if any
    const bool hasIdem = !idems.empty();
    std::vector<bool> can_be_first(d_table.order(), !hasIdem);
    for (const auto i : idems)
        can_be_first[i] = true;

    // count f(r,y)=r
    std::vector<size_t> rs(d_table.order(), 0);
    size_t mxrs = 0;
    for (auto row = d_table.order(); row--;) {
        for (auto col = d_table.order(); col--;)
            if (d_table.get(row, col) == row)
                rs[row]++;
        if (rs[row] > mxrs)
            mxrs = rs[row];
    }
    if (mxrs > 0) {
        for (auto row = d_table.order(); row--;)
            if (rs[row] != mxrs)
                can_be_first[row] = false;
    }
    size_t count_can_be_first = 0;
    size_t maxRow = -1;
    for (auto row = d_table.order(); row--;)
        if (!can_be_first[row])
            d_sat->addClause(~d_encoding->perm(row, 0));
        else {
            count_can_be_first++;
            maxRow = row;
        }

    if (count_can_be_first == 1) {
        d_statistics.unique1stRow->inc();
        d_fixed[0] = maxRow;
    }

    assert(count_can_be_first);
}

bool LexminSolver::test_sat_inc() {
    Minisat::vec<Minisat::Lit> assumps(1, mkLit(d_sat->fresh()));
    const auto &selector = assumps[0];
    d_encoding->encode_pos(d_assignments.back(), selector);
    const auto res = d_sat->solve(assumps);
    d_sat->addClause(res ? selector : ~selector);
    return res;
}

bool LexminSolver::test_sat_noinc() {
    d_sat = std::make_unique<SATSPC::MiniSatExt>();
    make_encoding();
    d_encoding->encode(d_assignments);
    const auto res = d_sat->solve();
    d_encoding.reset();
    d_sat.reset();
    return res;
}

void LexminSolver::calculate_budgets() {
    const auto n = d_table.order();
    size_t max_occurrences = 0; // max occurrences for non-fixed elements
    std::vector<size_t> max_occurrences_fixed(n, 0); // max occs. fixed elems
    std::vector<size_t> occurrences; // occurrences in current row
    for (auto row = n; row--;) {
        if (d_used[row]) // skip used rows
            continue;
        // count occurrences for each value
        occurrences.clear();
        occurrences.resize(n, 0);
        for (auto col = n; col--;)
            occurrences[d_table.get(row, col)]++;
        // update occurrence maxima
        for (auto value = n; value--;) {
            auto &mx = is_fixed(value) ? max_occurrences_fixed[value]
                                       : max_occurrences;
            mx = std::max(mx, occurrences[value]);
        }
    }
    d_output.comment(2) << "max occ. per row: " << max_occurrences << std::endl;
    d_rowBudget.clear();
    d_rowBudget.resize(n, max_occurrences); // defaults to max_occurrences
    // set up budgets for fixed values' images
    for (auto original_value = n; original_value--;)
        if (is_fixed(original_value)) {
            const auto image = d_fixed[original_value];
            d_rowBudget[image] = max_occurrences_fixed[original_value];
            d_output.comment(2) << "max occ. per row for " << image << ": "
                                << d_rowBudget[image] << std::endl;
        }
}

void LexminSolver::make_solution() {
    const auto n = d_table.order();
    d_solution = std::make_unique<BinaryFunction>(n);
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

