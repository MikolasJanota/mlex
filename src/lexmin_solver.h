/*
 * File:  lexmin_solver.h
 * Author:  mikolas
 * Created on:  Wed Dec 14 13:29:46 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include "auxiliary.h"       // for SATSPC
#include "binary_function.h" // for BinaryFunction
#include "comp_function.h"
#include "encoding.h"
#include "fixed_elem_map.h"
#include "invariants.h"
#include "options.h"
#include <cstddef> // for size_t
#include <iosfwd>  // for ostream
#include <memory>
#include <optional> // for optional
#include <set>      // for set
#include <utility>  // for pair
#include <vector>

class RowBudgets;
class IBudget;
class StatisticsManager;

class LexminSolver {
  public:
    LexminSolver(Output &output, const BinaryFunction &table);
    virtual ~LexminSolver();

    /* Run the solver */
    void solve();

    /* Make solution as compact BinaryFunction */
    BinaryFunction *make_solution();

    /* Make solution as compact fun representation */
    CompFunction make_solution_comp();

    /* Only if d_options.diagonal */
    void set_diag(const std::vector<size_t> &diag) {
        assert(diag.size() == d_table.order());
        d_diag = diag;
    }

    /* Only if d_options.diagonal */
    void set_diag_(std::vector<size_t> &diag) {
        assert(diag.size() == d_table.order());
        d_diag = std::move(diag);
    }

  private:
    Output &d_output;
    const Options &d_options;
    StatisticsManager &d_statistics;
    const BinaryFunction &d_table;

    bool d_is_solved = false;
    std::vector<size_t> d_diag;

    std::unique_ptr<SATSPC::MiniSatExt> d_sat;
    std::unique_ptr<Encoding> d_encoding;

    std::vector<Encoding::Assignment> d_assignments;

    FixedElemMap d_fixed;
    std::vector<bool> d_used;

    /** cells already fixed in the output table **/
    std::unique_ptr<BinaryFunction> d_fixed_cells;

    Invariants d_invariants;
    std::unique_ptr<RowBudgets> d_budgets;

    inline std::ostream &comment(int level = 0) {
        return d_output.comment(level);
    }

    inline std::ostream &ccomment(int level = 0) {
        return d_output.ccomment(level);
    }

    void make_encoding();
    bool test_sat(const Encoding::Assignment &asg);
    bool test_sat(const std::pair<size_t, size_t> &cell,
                  const std::vector<size_t> &vals);
    bool test_sat_inc(const Encoding::Assignment &asg);

    /* try to infer additional constraints on the first row */
    void opt1stRow();
    void calculate_budgets_row_tot();
    void calculate_budgets_col();
    void run_diagonal();
    size_t id_row_elements(size_t dst_row, size_t src_row);
    void closure_fixed();
    void calculate_diagonal(size_t max_idem_reps,
                            const DiagInvariants &di_orig);
    void mark_used_rows(const Invariants::Info &rows, size_t current_row);

    // returns whether updates should be updated
    bool process_invariant(const InvariantVector &invv, size_t current_row);
    size_t find_value(Encoding::Assignment &a, IBudget &budget,
                      const std::optional<size_t> &last_val);
    size_t find_value_unsat_sat(Encoding::Assignment &a, IBudget &budget,
                                const std::optional<size_t> &last_val);
    size_t find_value_sat_unsat(Encoding::Assignment &a, IBudget &budget,
                                const std::optional<size_t> &last_val);
    size_t find_value_bin(Encoding::Assignment &a, IBudget &budget,
                          const std::optional<size_t> &last_val);
    size_t find_value_bin2(Encoding::Assignment &a, IBudget &budget,
                           const std::optional<size_t> &last_val);

    std::vector<size_t> d_last_permutation;
    std::vector<size_t> d_inv_last_permutation;
    void make_last_permutation();
    size_t get_preimage(size_t i) const;
    size_t get_image(size_t i) const;
    size_t get_val(size_t row, size_t col) const;
    std::ostream &show_permutation(std::ostream &out);
    std::ostream &show_diag(std::ostream &out);
};
