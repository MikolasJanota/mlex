/*
 * File:  lexmin_solver.h
 * Author:  mikolas
 * Created on:  Wed Dec 14 13:29:46 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include "encoding.h"
#include "invariants.h"
#include "minisat_ext.h"
#include "options.h"
#include <memory>
#include <vector>

class LexminSolver {
  public:
    LexminSolver(Output &output, const BinaryFunction &table)
        : d_output(output), d_options(output.d_options),
          d_statistics(output.d_statistics), d_table(table),
          d_invariants(output, table) {}
    void solve();

    void print_solution(std::ostream &output);

    std::unique_ptr<BinaryFunction> &solution() { return d_solution; }

  private:
    Output &d_output;
    const Options &d_options;
    StatisticsManager &d_statistics;
    const BinaryFunction &d_table;
    std::unique_ptr<BinaryFunction> d_solution;
    std::unique_ptr<SATSPC::MiniSatExt> d_sat;
    std::unique_ptr<Encoding> d_encoding;

    std::vector<Encoding::Assignment> d_assignments;
    std::vector<size_t> d_row_budget;
    std::vector<size_t> d_col_budget;
    std::vector<size_t> d_total_budget;

    std::vector<size_t> d_fixed;
    std::vector<bool> d_used;
    bool is_fixed(size_t i) const { return d_fixed[i] < d_table.order(); }
    std::optional<size_t> d_0preimage;

    Invariants d_invariants;

    inline std::ostream &comment(int level = 0) {
        return d_output.comment(level);
    }

    inline std::ostream &ccomment(int level = 0) {
        return d_output.ccomment(level);
    }

    void make_encoding();
    bool test_sat();
    bool test_sat_noinc();
    bool test_sat_inc();

    /*  try to infer additional constraints on the first row */
    void opt1stRow();
    void make_solution();
    void calculate_budgetsRowTot();
    void calculate_budgetsCol();
    void mark_used_rows(const Invariants::Info &rows, size_t current_row);

    // returns whether updates should be updated
    bool process_invariant(const InvariantVector &invv, size_t current_row);
};
