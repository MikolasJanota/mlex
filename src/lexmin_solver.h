/*
 * File:  lexmin_solver.h
 * Author:  mikolas
 * Created on:  Wed Dec 14 13:29:46 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include "encoding.h"
#include "minisat_ext.h"
#include "options.h"
#include <memory>
#include <vector>

class LexminSolver {
  public:
    LexminSolver(Output &output, const BinaryFunction &table)
        : d_output(output), d_options(output.d_options),
          d_statistics(output.d_statistics), d_table(table) {}
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
    std::vector<size_t> d_rowBudget;

    bool test_sat();
    bool test_sat_noinc();
    bool test_sat_inc();

    void make_solution();
    void calculate_basic_budgets();
};
