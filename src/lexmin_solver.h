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
    LexminSolver(const Options &options, const BinaryFunction &table)
        : d_options(options), d_table(table) {}
    void solve();

    void print_solution(std::ostream &output);

    std::unique_ptr<BinaryFunction> &solution() { return d_solution; }

  private:
    const Options &                     d_options;
    double                              d_total_sat_time = 0;
    const BinaryFunction &              d_table;
    std::unique_ptr<BinaryFunction>     d_solution;
    std::unique_ptr<SATSPC::MiniSatExt> d_sat;
    std::unique_ptr<Encoding>           d_encoding;
    bool                                test_sat(const Encoding::Assignment &);
    bool test_sat_noinc(const Encoding::Assignment &);
    bool test_sat_inc(const Encoding::Assignment &);
    std::vector<Encoding::Assignment> d_assignments;

    void make_solution();
    void print_gap(std::ostream &output);
    void print_mace(std::ostream &output);

    std::ostream &comment() {
        return std::cout << (d_options.mace_format ? '%' : '#') << " ";
    }
};
