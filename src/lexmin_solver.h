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

class Budgets;

class LexminSolver {
  public:
    LexminSolver(Output &output, const BinaryFunction &table);
    virtual ~LexminSolver();
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

    std::vector<size_t> d_fixed;
    std::vector<bool> d_used;
    bool is_fixed(size_t i) const { return d_fixed[i] < d_table.order(); }
    std::optional<size_t> d_0preimage;

    Invariants d_invariants;
    std::unique_ptr<Budgets> d_budgets;

    inline std::ostream &comment(int level = 0) {
        return d_output.comment(level);
    }

    inline std::ostream &ccomment(int level = 0) {
        return d_output.ccomment(level);
    }

    void make_encoding();
    bool test_sat();
    bool test_sat(const std::pair<size_t, size_t> &cell,
                  const std::vector<size_t> &vals);
    bool test_sat_noinc();
    bool test_sat_inc();

    /*  try to infer additional constraints on the first row */
    void opt1stRow();
    void make_solution();
    void calculate_budgets_row_tot();
    void calculate_budgets_col();
    void mark_used_rows(const Invariants::Info &rows, size_t current_row);

    // returns whether updates should be updated
    bool process_invariant(const InvariantVector &invv, size_t current_row);
    size_t find_value(std::optional<size_t> last_val);
    size_t find_value_unsat_sat(std::optional<size_t> last_val);
    size_t find_value_sat_unsat(std::optional<size_t> last_val);

    std::vector<size_t> d_last_permutation;
    std::vector<size_t> d_inv_last_permutation;
    void make_last_permutation();
    size_t get_preimage(size_t i) const;
    size_t get_image(size_t i) const;
    size_t get_val(size_t row, size_t col) const;
};
