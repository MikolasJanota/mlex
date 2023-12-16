/*
 * File:  lexmin_solver.cpp
 * Author:  mikolas
 * Created on:  Wed Dec 14 13:29:54 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "lexmin_solver.h"
#include "auxiliary.h"
#include "encoding.h"
#include "immutable_vector.h"
#include "invariants.h"
#include "minisat/core/SolverTypes.h"
#include "minisat_ext.h"
#include "options.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <tuple> // for get
#include <vector>
using SATSPC::Lit;
using SATSPC::mkLit;
using std::endl;
/* #define SOLVER_TRACING */
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

class IBudget {
  public:
    virtual bool has_budget(size_t val) const = 0;
};

class TrivBudget : public IBudget {
  public:
    virtual bool has_budget(size_t) const { return true; }
};

class DiagBudget : public IBudget {
  public:
    DiagBudget(size_t ord, size_t reps, size_t idem_count)
        : d_budgets(ord, reps), d_idem_count(idem_count) {}
    virtual ~DiagBudget() {}

    virtual bool has_budget(size_t val) const {
        return (val != d_row || d_idem_count > 0) && d_budgets[val] > 0;
    }

    void dec(size_t val) {
        assert(d_budgets[val] > 0);
        d_budgets[val]--;
        if (val == d_row) {
            assert(d_idem_count > 0);
            d_idem_count--;
        }
    }

    void set_row(size_t row) { d_row = row; }

    std::vector<size_t> d_budgets;
    size_t d_idem_count;

  private:
    size_t d_row;
};

#define CHECKED_DECREASE(v)                                                    \
    do {                                                                       \
        assert((v) > 0);                                                       \
        (v)--;                                                                 \
    } while (false);

class RowBudgets {
  public:
    RowBudgets(const BinaryFunction &table) : d_table(table) {}

    bool has_budget(size_t col, size_t val) const {
        assert(col < d_table.order());
        assert(val < d_table.order());
        return d_current_row_budget[val] > 0 && d_total_budget[val] > 0 &&
               d_cols_budget[col][val] > 0;
    }

    void dec_budget(size_t col, size_t val) {
        assert(col < d_table.order());
        assert(val < d_table.order());

        CHECKED_DECREASE(d_current_row_budget[val]);
        CHECKED_DECREASE(d_total_budget[val]);
        CHECKED_DECREASE(d_cols_budget[col][val]);
    }

    std::ostream &print(std::ostream &out, size_t col, size_t val) const {
        return out << "(budg r/c/t: " << d_current_row_budget[val] << " "
                   << d_cols_budget[col][val] << " " << d_total_budget[val]
                   << ")";
    }

    void reset_cur_row_budget() {
        assert(d_all_budget.size() == d_table.order());
        d_current_row_budget = d_all_budget;
    }

    void reset_cols_budget() {
        assert(d_col_budget.size() == d_table.order());
        d_cols_budget.assign(d_table.order(), d_col_budget);
    }

    std::vector<std::vector<size_t>> d_cols_budget;

    /* Currently we have 3 types of budgets for rows. Rows that contain an
     * idempotent, rows that do not contain idempotent and all rows.
     * The idea is that once we find out  whether the current row has or does
     * not have an idempotent, we can refine the budget. */
    typedef std::vector<size_t> RowBudget;
    RowBudget d_idem_budget;
    RowBudget d_noidem_budget;
    RowBudget d_all_budget;

    void refine_budget(bool idem) {
        const RowBudget &newb = idem ? d_idem_budget : d_noidem_budget;
        for (size_t i = d_current_row_budget.size(); i--;) {
            const auto d = d_all_budget[i] - newb[i];
            assert(d_all_budget[i] >= newb[i]);
            assert(d_current_row_budget[i] >= d);
            d_current_row_budget[i] -= d;
        }
    }

    std::vector<size_t> d_col_budget;
    std::vector<size_t> d_current_row_budget;
    std::vector<size_t> d_total_budget;

  private:
    const BinaryFunction &d_table;
};

class RowBudgeter : public IBudget {
  public:
    RowBudgeter(RowBudgets *paren, bool trivial)
        : d_paren(paren), d_triv(trivial), d_col(-1) {}

    void set_col(size_t col) { d_col = col; }

    virtual bool has_budget(size_t val) const override {
        return d_triv || d_paren->has_budget(d_col, val);
    }

  private:
    RowBudgets *d_paren;
    bool d_triv;
    size_t d_col;
};

LexminSolver::LexminSolver(Output &output, const BinaryFunction &table)
    : LexminSolverBase(output, table), d_fixed(output, table.order()),
      d_invariants(output, table) {}

LexminSolver::~LexminSolver() {}

size_t LexminSolver::find_value(Encoding::Assignment &asg, IBudget &budget,
                                const std::optional<size_t> &last_val) {
    /* const auto n = d_table.order(); */
    auto &[row, col, val] = asg;
    if (d_fixed_cells) {
        if (d_fixed_cells->is_set(row, col)) {
            val = d_fixed_cells->get(row, col);
            TRACE(comment(3)
                      << "fixed cell: (" << row << "," << col << ")=" << val;);
            return val;
        }
    }

    /* ? (row == 0 && col < n / 2 ? SearchType::lin_us : SearchType::bin2) */
    const auto st = d_options.search_type == SearchType::adaptive
                        ? (row == 0 ? SearchType::lin_us : SearchType::bin2)
                        : d_options.search_type;

    switch (st) {
    case SearchType::lin_us: return find_value_unsat_sat(asg, budget, last_val);
    case SearchType::lin_su: return find_value_sat_unsat(asg, budget, last_val);
    case SearchType::bin: return find_value_bin(asg, budget, last_val);
    case SearchType::bin2: return find_value_bin2(asg, budget, last_val);
    case SearchType::adaptive:
        assert(0);
        return find_value_bin2(asg, budget, last_val);
    }
    assert(false);
    return -1;
}

size_t LexminSolver::find_value_bin2(Encoding::Assignment &asg, IBudget &budget,
                                     const std::optional<size_t> &last_val) {
    const auto n = d_table.order();
    auto &[row, col, cur_val] = asg;
    const auto cell = std::make_pair<>(row, col);

    TRACE(comment(3) << "(" << row << " " << col << ") :";);
    size_t ub = last_val ? *last_val : n; // upper bound
    TRACE(comment(3) << "(iub:" << ub << ") ";);
    std::vector<size_t> a, b;
    std::vector<size_t> *vals = &a, *top = &b;

    for (size_t v = 0; v < ub; v++) {
        if (!d_options.budgeting || budget.has_budget(v))
            vals->push_back(v);
    }

    TRACE(print_set(comment(4) << "vals ", *vals););
    while (!vals->empty()) {
        assert(top->empty());

        if (!test_sat(cell, *vals)) {
            break;
        }
        ub = std::min(ub, get_val(row, col));
        TRACE(comment(3) << "(ub:" << ub << ") ";);

        const auto split = vals->size() / 2 + vals->size() % 2;
        for (size_t h = split; h < vals->size(); h++) {
            const auto val = vals->at(h);
            if (val < ub)
                top->push_back(val);
        }
        vals->resize(split);
        while (!vals->empty() && vals->back() >= ub)
            vals->pop_back();

        if (top->empty())
            continue; // ub fell into first half

        if (test_sat(cell, *vals)) {
            assert(!d_last_permutation.empty());
            ub = std::min(ub, get_val(row, col));
            TRACE(comment(3) << "(ub:" << ub << ") ";);
            while (!vals->empty() && vals->back() >= ub)
                vals->pop_back();
        } else {
            std::swap(vals, top);
        }
        top->clear();
    }

    cur_val = ub;
    TRACE(comment(3) << "val:" << cur_val;);
    d_encoding->encode_pos(asg, SATSPC::lit_Undef);
    return cur_val;
}

size_t LexminSolver::find_value_bin(Encoding::Assignment &asg, IBudget &budget,
                                    const std::optional<size_t> &last_val) {
    const auto n = d_table.order();
    auto &[row, col, cur_val] = asg;
    const auto cell = std::make_pair<>(row, col);

    TRACE(comment(3) << "(" << row << " " << col << ") :";);
    size_t ub = last_val ? *last_val : n; // upper bound
    TRACE(comment(3) << "(iub:" << ub << ") ";);
    std::vector<size_t> a, b;
    std::vector<size_t> *vals = &a, *top = &b;

    for (size_t v = 0; v < ub; v++)
        if (budget.has_budget(v))
            vals->push_back(v);

    while (!vals->empty()) {
        assert(top->empty());
        TRACE(print_set(comment(4) << "vals ", *vals););

        const auto split = vals->size() / 2 + vals->size() % 2;
        for (size_t h = split; h < vals->size(); h++)
            top->push_back(vals->at(h));
        vals->resize(split);

        if (test_sat(cell, *vals)) {
            assert(!d_last_permutation.empty());
            ub = std::min(ub, get_val(row, col));
            TRACE(comment(3) << "(ub:" << ub << ") ";);
            while (!vals->empty() && vals->back() >= ub)
                vals->pop_back();
        } else {
            std::swap(vals, top);
        }
        top->clear();
    }

    cur_val = ub;
    TRACE(comment(4) << "val: " << cur_val;);
    d_encoding->encode_pos(asg, SATSPC::lit_Undef);
    return cur_val;
}

size_t
LexminSolver::find_value_sat_unsat(Encoding::Assignment &asg, IBudget &budget,
                                   const std::optional<size_t> &last_val) {

    const auto n = d_table.order();
    auto &[row, col, cur_val] = asg;
    const auto cell = std::make_pair<>(row, col);

    TRACE(comment(3) << "(" << row << " " << col << ") :";);
    size_t ub = last_val ? *last_val : n; // upper bound
    TRACE(comment(3) << "(iub:" << ub << ") ";);
    std::vector<size_t> vals;
    for (size_t v = 0; v < ub; v++) {
        if (budget.has_budget(v))
            vals.push_back(v);
    }
    while (!vals.empty() && test_sat(cell, vals)) {
        assert(!d_last_permutation.empty());
        ub = std::min(ub, get_val(row, col));
        TRACE(comment(3) << "(ub:" << ub << ") ";);
        while (!vals.empty() && vals.back() >= ub)
            vals.pop_back();
    }
    cur_val = ub;
    TRACE(comment(4) << "val: " << cur_val;);
    d_encoding->encode_pos(asg, SATSPC::lit_Undef);
    return ub;
}

size_t
LexminSolver::find_value_unsat_sat(Encoding::Assignment &asg, IBudget &budget,
                                   const std::optional<size_t> &last_val) {
    auto &[row, col, cur_val] = asg;
    TRACE(comment(3) << "(" << row << " " << col << ") :";);

    for (;; cur_val++) {
        assert(cur_val < d_table.order());
        if (!budget.has_budget(cur_val))
            continue;
        const bool free = last_val && *last_val == cur_val;
        if (free) {
            TRACE(ccomment(3) << " " << cur_val << ":FREE";);
            d_encoding->encode_pos(asg, SATSPC::lit_Undef);
        }
        if (free || test_sat(asg))
            return cur_val;
    }
}

/* The output table will contain max_idem_reps 0's at the beginning. If
 * max_idem_reps > 0 these are the only 0's.*/
void LexminSolver::calculate_diagonal(size_t max_idem_reps,
                                      const DiagInvariants &di_orig) {
    using std::max;
    assert(d_fixed.fixed_count() <= 1); // fixed first idem
    const auto n = d_table.order();

    // statistics for non-fixed
    size_t max_reps = 0;
    size_t idem_count = 0;
    for (size_t i = n; i--;) {
        if (d_fixed.is_fixed_src(i))
            continue;
        max_reps = max(max_reps, di_orig.get_reps(i));
        if (d_table.get(i, i) == i)
            idem_count++;
    }

    // set up budgeting for the diagonal
    DiagBudget budg(n, max_reps, idem_count);
    if (max_idem_reps > 0)
        budg.d_budgets[0] = 0;

    // place max_idem_reps 0's at beginning
    for (auto i = max_idem_reps; i--;) {
        d_fixed_cells->set(i, i, 0);
        d_encoding->encode_pos({i, i, 0}, SATSPC::lit_Undef);
    }

    // calculate the rest of the diagonal
    for (size_t i = max_idem_reps; i < n; i++) {
        Encoding::Assignment asg{i, i, 0};
        const auto last_val = d_last_permutation.empty()
                                  ? std::nullopt
                                  : std::make_optional(get_val(i, i));
        budg.set_row(i);
        const auto new_val = find_value(asg, budg, last_val);
        d_fixed_cells->set(i, i, new_val);
        budg.dec(new_val);
        TRACE(ccomment(3) << std::endl;);
    }
    d_statistics.satCalls->print(comment(2) << "diag:") << std::endl;
}

void LexminSolver::run_diagonal() {
    using std::max;
    const auto n = d_table.order();

    // process invariants in the input table
    DiagInvariants di_orig(d_output, n);
    for (size_t i = n; i--;) {
        di_orig.add(i);
        di_orig.set(i, d_table.get(i, i));
    }
    di_orig.calculate();
    size_t max_idem_reps = 0;  // max reps of an idem
    std::vector<size_t> idems; // list of idems
    for (size_t i = n; i--;) {
        if (di_orig.get_loop(i) == 1) {
            // idempotent
            max_idem_reps = max(max_idem_reps, di_orig.get_reps(i));
            idems.push_back(i);
        }
    }

    // try identifying candidates for 0
    if (!idems.empty()) {
        // idempotents that appear the most times
        std::set<size_t> candidates;
        for (auto i : idems)
            if (di_orig.get_reps(i) == max_idem_reps)
                candidates.insert(i);
        if (candidates.size() == 1) { // top left corner fixed
            d_statistics.uniqueDiag1Elem->inc();
            d_fixed.set(first(candidates), 0);
            comment(2) << first(candidates) << " fixed to 0 (diagonal)"
                       << std::endl;
        }
        // prohibit non-candidates to be mapped to 0
        for (auto i = n; i--;)
            if (!contains(candidates, i))
                d_sat->addClause(~d_encoding->perm(i, 0));
    }

    if (d_diag.empty()) {
        calculate_diagonal(max_idem_reps, di_orig);
    } else {
        for (size_t i = 0; i < n; ++i) {
            const auto val = d_diag[i];
            comment(3) << i << ":" << val << " by diag file" << std::endl;
            assert(i >= max_idem_reps || val == 0);
            d_fixed_cells->set(i, i, val);
            d_encoding->encode_pos({i, i, val}, SATSPC::lit_Undef);
        }
    }
    show_diag(comment(2) << "diag:") << std::endl;

    // calculate invariants for the destination table's diagonal
    DiagInvariants di_new(d_output, n);
    for (size_t i = 0; i < n; i++) {
        di_new.add(i);
        di_new.set(i, d_fixed_cells->get(i, i));
    }
    di_new.calculate();
    di_new.calc_inverse();

    // identify elements based on corresponding invariants
    for (size_t src = 0; src < n; src++) {
        const auto &dest_corresp = di_new.get_elems(di_orig.get_invariant(src));
        // disable permuting to nonmatching elements
        for (size_t dest = n; dest--;)
            if (!contains(dest_corresp, dest))
                d_sat->addClause(~d_encoding->perm(src, dest));
        // handle the singleton case
        if (dest_corresp.size() == 1) {
            d_statistics.uniqueDiagElem->inc();
            if (d_fixed.set(src, first(dest_corresp))) {
                comment(2) << src << " fixed to " << d_fixed.src2dst(src)
                           << " (diag)" << std::endl;
            }
        }
    }

    closure_fixed();
}

void LexminSolver::closure_fixed() {
    for (auto row : d_fixed.fixed_src_elements())
        for (auto col : d_fixed.fixed_src_elements()) {
            const auto src_val = d_table.get(row, col);
            if (!d_fixed.is_fixed_src(src_val)) {
                continue;
            }
            const auto rdst = d_fixed.src2dst(row);
            const auto cdst = d_fixed.src2dst(col);
            const auto vdst = d_fixed.src2dst(src_val);
            if (d_fixed_cells->is_set(rdst, cdst)) {
                assert(d_fixed_cells->get(rdst, cdst) == vdst);
                continue;
            }
            d_fixed_cells->set(rdst, cdst, vdst);
            d_statistics.inferredCells->inc();
            comment(2) << "fixing (" << rdst << "," << cdst << ")=" << vdst
                       << std::endl;
        }
}

void LexminSolver::solve() {
    const auto n = d_table.order();
    const auto budgeting = d_options.budgeting;
    d_used.resize(n, false);

    if (d_options.diagonal || d_options.opt1stRow || d_options.id_elements)
        d_fixed_cells = std::make_unique<BinaryFunction>(n);

    if (d_options.invariants || d_options.inv_ord)
        d_invariants.calculate();

    std::unique_ptr<ColorInvariantCalculator> ccalc;

    if (d_options.color) {
        d_colors = std::make_unique<ColorInvariantManager>(d_output, d_table);
    }

    d_sat = std::make_unique<SATSPC::MiniSatExt>();
    d_encoding = std::make_unique<Encoding>(d_output, *d_sat, d_table);
    if (d_options.opt1stRow && !d_options.diagonal)
        opt1stRow();
    if (d_options.inv_ord)
        enc_inv_ord();
    d_encoding->encode_bij();

    if (d_options.diagonal)
        run_diagonal();

    if (budgeting) {
        d_budgets = std::make_unique<RowBudgets>(d_table);
        calculate_budgets_row_tot();
        calculate_budgets_col();
        d_budgets->reset_cols_budget();
    }

    d_assignments.clear();
    d_assignments.reserve(n * n);

    std::vector<size_t> row_vals(n, -1);
    InvariantCalculator calc(n);

    const auto start_time = read_cpu_time();
    for (size_t row = 0; row < n; row++) {
        comment(2) << "row:" << row << " "
                   << SHOW_TIME(read_cpu_time() - start_time) << std::endl;

        const auto fixed_before = d_fixed.fixed_count();

        if (d_options.simp_sat_row)
            d_sat->simplify();

        if (budgeting) {
            d_budgets->reset_cur_row_budget();
            if (d_options.diagonal && d_options.budget_idem)
                d_budgets->refine_budget(d_fixed_cells->get(row, row) == row);
        }

        if (d_options.invariants) {
            calc.set_row(row);
        }

        RowBudgeter budgeter(d_budgets.get(), !d_budgets);
        for (size_t col = 0; col < n; col++) {
            budgeter.set_col(col);
            d_assignments.push_back({row, col, 0});
            const auto cur_val =
                find_value(d_assignments.back(), budgeter,
                           d_last_permutation.empty()
                               ? std::nullopt
                               : std::make_optional(get_val(row, col)));
            assert(cur_val == std::get<2>(d_assignments.back()));
            row_vals[col] = cur_val;

            if (d_budgets) {
                d_budgets->dec_budget(col, cur_val);
                if (!d_options.diagonal && d_options.budget_idem && col == row)
                    d_budgets->refine_budget(cur_val == row);
            }

            if (d_options.invariants)
                calc.set_val(col, cur_val);

            if (d_fixed_cells)
                d_fixed_cells->set(row, col, cur_val);

            TRACE(ccomment(3) << std::endl;);
        }

        bool update_budgets = false;

        if (d_options.invariants && row + 1 < n) {
            Looping lc(d_output, row_vals);
            Distances dc(d_output, row_vals, row);
            for (auto col = n; col--;) {
                calc.add_loop(lc.calc_loop(col));
                if (d_output.d_options.distance_invariant)
                    calc.add_distance(dc.calc_distance(col));
            }
            update_budgets = process_invariant(calc.make_ivec(), row);
        }

        if (row == 0 && d_fixed.is_fixed_dst(0)) {
            d_used[d_fixed.dst2src(0)] = true; // mark preimage of row 0 as used
            update_budgets = true;
        }

        if (row + 1 < n) {
            if (d_fixed.is_fixed_dst(row))
                id_row_elements(row, d_fixed.dst2src(row));
            if (d_options.color)
                process_invariant_color(row);
            if (fixed_before != d_fixed.fixed_count() && d_fixed_cells)
                closure_fixed();
            if (update_budgets && budgeting)
                calculate_budgets_row_tot();
        }
    }
    d_is_solved = true;
    if (!d_last_permutation.empty())
        show_permutation(comment(2) << "perm:", d_last_permutation)
            << std::endl;
}

size_t LexminSolver::id_row_elements(size_t dst_row, size_t src_row) {
    if (!d_options.id_elements)
        return 0;
    comment(2) << "id_elements:" << dst_row << "<-" << src_row << std::endl;
    assert(d_fixed.src2dst(src_row) == dst_row);
    const auto n = d_table.order();
    // calculate invariants in the source table's row
    DiagInvariants inv_src(d_output, n);
    for (auto col = n; col--;) {
        if (!d_fixed.is_fixed_src(col))
            inv_src.add(col);
        inv_src.set(col, d_table.get(src_row, col));
    }

    if (inv_src.elems().empty())
        return 0; // everything already fixed
    inv_src.calculate();

    // calculate invariants in the destination table's row
    DiagInvariants inv_dst(d_output, n);
    for (auto col = n; col--;) {
        if (!d_fixed.is_fixed_dst(col))
            inv_dst.add(col);
        inv_dst.set(col, d_fixed_cells->get(dst_row, col));
    }
    inv_dst.calculate();
    inv_dst.calc_inverse();

    size_t rv = 0;
    // identify elements based on corresponding invariants
    for (size_t src_elem : inv_src.elems()) {
        // get elements in dst with the same invariant as src_elem
        const auto &dst_corresp =
            inv_dst.get_elems(inv_src.get_invariant(src_elem));
        // disable permuting to nonmatching elements
        for (size_t dst_elem = n; dst_elem--;)
            if (!contains(dst_corresp, dst_elem))
                d_sat->addClause(~d_encoding->perm(src_elem, dst_elem));
        // handle the singleton case
        if (dst_corresp.size() == 1) {
            rv++;
            if (d_fixed.set(src_elem, first(dst_corresp))) {
                d_statistics.uniqueRowElem->inc();
                comment(2) << src_elem << " fixed to "
                           << d_fixed.src2dst(src_elem) << " (idrow)"
                           << std::endl;
            }
        }
    }
    return rv;
}

std::ostream &LexminSolver::show_diag(std::ostream &out) {
    const auto n = d_table.order();
    out << '[';
    for (size_t i = 0; i < n; i++)
        out << (i ? ", " : " ") << d_fixed_cells->get(i, i);
    return out << ']';
}

void LexminSolver::make_last_permutation() {
    const auto n = d_table.order();
    const auto first = d_last_permutation.empty();
    if (first) {
        d_last_permutation.assign(n, -1);
        d_inv_last_permutation.assign(n, -1);
    }
    assert(d_last_permutation.size() == n);
    assert(d_inv_last_permutation.size() == n);
    for (auto dom = n; dom--;) {
        bool found =
            !first &&
            d_sat->eval_lit(d_encoding->perm(dom, d_last_permutation[dom])) ==
                SATSPC::l_True;
        for (auto rng = n; !found && rng--;)
            if (d_sat->eval_lit(d_encoding->perm(dom, rng)) == SATSPC::l_True) {
                d_last_permutation[dom] = rng;
                /* std::cerr << dom << "->" << rng << std::endl; */
                d_inv_last_permutation[rng] = dom;
                found = true;
            }
        assert(found);
    }
}

/* size_t LexminSolver::get_preimage(size_t i) { */
/*     const auto n = d_table.order(); */
/*     for (auto j = n; j--;) */
/*         if (d_sat->eval_lit(d_encoding->perm(j, i)) == SATSPC::l_True) */
/*             return j; */
/*     for (size_t r = 0; r < n; r++) { */
/*         comment(2) << r; */
/*         for (size_t s = 0; s < n; s++) { */
/*             const auto v = d_sat->eval_lit(d_encoding->perm(r, s)); */
/*             ccomment(2) << " " */
/*                         << (v == SATSPC::l_True */
/*                                 ? "+" */
/*                                 : (v == SATSPC::l_False ? "-" : "?")); */
/*         } */
/*         ccomment(2) << std::endl; */
/*     } */
/*     assert(false); */
/*     return -1; */
/* } */

size_t LexminSolver::get_val(size_t row, size_t col) const {
    const auto row_preim = d_inv_last_permutation[row];
    const auto col_preim = d_inv_last_permutation[col];
    const auto val = d_table.get(row_preim, col_preim);
    return d_last_permutation[val];
}

void LexminSolver::match_lines_color(size_t src_row, size_t dst_row) {
    const auto n = d_table.order();
    d_output.comment(2) << "[col] match lines:" << src_row << "->" << dst_row
                        << endl;
    if (d_fixed.set(src_row, dst_row)) {
        comment(2) << "[col] " << src_row << " fixed to " << dst_row << endl;
        d_statistics.uniqueColInv->inc();
    }
    ColorInvariantCalculator srcc(d_output, d_colors->d_color_count,
                                  d_colors->d_colors_src);
    srcc.set_row(src_row);
    for (size_t col = 0; col < n; ++col)
        srcc.set_val(col, d_table.get(src_row, col));
    ColorInvariantCalculator::InvMap src_map;
    srcc.make_map(src_map);

    ColorInvariantCalculator dstc(d_output, d_colors->d_color_count,
                                  d_colors->d_colors_dst);
    dstc.set_row(dst_row);
    for (size_t col = 0; col < n; ++col)
        dstc.set_val(col, d_fixed_cells->get(dst_row, col));
    ColorInvariantCalculator::InvMap dst_map;
    dstc.make_map(dst_map);
    size_t dcount = 0;

    for (const auto &[inv, dst_elems] : dst_map) {
        const auto &src_elems = src_map.at(inv);
        print_vec(d_output.comment(3) << "inv:", inv) << std::endl;
        print_vec(d_output.comment(3) << "   dst:", dst_elems) << std::endl;
        print_vec(d_output.comment(3) << "   src:", src_elems) << std::endl;
        assert(src_elems.size() == dst_elems.size());
        for (const auto dst_elem : dst_elems)
            for (size_t src_elem = 0; src_elem < n; src_elem++)
                if (!contains(src_elems, src_elem))
                    if (try_disallow(src_elem, dst_elem))
                        dcount++;
        for (const auto src_elem : src_elems)
            for (size_t dst_elem = 0; dst_elem < n; dst_elem++)
                if (!contains(dst_elems, dst_elem))
                    if (try_disallow(src_elem, dst_elem))
                        dcount++;

        if (src_elems.size() == 1) {
            const auto se = *(src_elems.begin());
            const auto de = *(dst_elems.begin());
            if (d_fixed.set(se, de)) {
                comment(2) << "[col] " << se << " fixed to " << de << endl;
                d_statistics.uniqueColInv->inc();
            }
        }
    }
}

bool LexminSolver::process_invariant_color(size_t current_row) {
    size_t old = -1;
    bool changed = false;
    do {
        old = d_colors->d_color_count;
        changed |= core_process_invariant_color(current_row);
    } while (0 && d_colors->d_color_count != old);
    return changed;
}

bool LexminSolver::core_process_invariant_color(size_t current_row) {
    const auto n = d_table.order();
    d_colors->add_row(current_row, *(d_fixed_cells.get()));
    if (static_cast<int>(d_colors->d_color_count) > d_statistics.maxCol->get())
        d_statistics.maxCol->set(d_colors->d_color_count);

    size_t dcount = 0;
    for (const auto &[inv, dst_rows] : d_colors->d_dst_row_color_invariants) {
        const auto &src_rows = d_colors->d_src_row_color_invariants.at(inv);
        assert(dst_rows.size() <= src_rows.size());
        d_output.comment(3) << "inv:" << std::endl;
        d_colors->print(inv);
        print_vec(d_output.comment(3) << "   dst:", dst_rows) << std::endl;
        print_vec(d_output.comment(3) << "   src:", src_rows) << std::endl;
        for (const auto dst_row : dst_rows)
            for (size_t src_row = 0; src_row < n; src_row++)
                if (!contains(src_rows, src_row))
                    if (try_disallow(src_row, dst_row))
                        dcount++;
        if (dst_rows.size() == src_rows.size()) {
            for (const auto src_row : src_rows)
                for (size_t dst_row = 0; dst_row < n; dst_row++)
                    if (!contains(dst_rows, dst_row))
                        if (try_disallow(src_row, dst_row))
                            dcount++;
            if (dst_rows.size() == 1)
                match_lines_color(*src_rows.begin(), *dst_rows.begin());
        }
    }

    return dcount > 0;
}

void LexminSolver::allowonly(size_t from, size_t to) {
    const auto n = d_table.order();
    for (size_t i = 0; i < n; i++)
        if (i != to)
            try_disallow(from, i);
    for (size_t i = 0; i < n; i++)
        if (i != from)
            try_disallow(i, to);
}

bool LexminSolver::try_disallow(size_t from, size_t to) {
    /* comment(4) << "try_disallow " << from << "->" << to << "@" */
    /*            << d_colors->d_allowed_mapping.allowed(from, to) << endl; */
    const auto rv = d_colors->d_allowed_mapping.disallow(from, to);
    if (rv) {
        d_sat->addClause(~d_encoding->perm(from, to));
        d_statistics.disallowed->inc();
        comment(4) << "disallowed " << from << "->" << to << endl;
    }
    return rv;
}
bool LexminSolver::process_invariant(const InvariantVector &invv,
                                     size_t current_row) {
    auto &info = d_invariants.get(invv);
    info.used++;
    const bool used_up = info.used == info.original_rows.size();
    if (used_up)
        mark_used_rows(info, current_row);
    return used_up;
}

void LexminSolver::mark_used_rows(const Invariants::Info &info,
                                  size_t current_row) {
    const auto &rows = info.original_rows;
    // marked as used for budgeting purposes
    for (auto k : rows)
        d_used[k] = true;

    TRACE(print_set(comment(3) << "used ", rows) << std::endl;);

    // disable further rows to be mapped to the used ones
    for (auto k : rows)
        for (size_t j = current_row + 1; j < d_table.order(); j++)
            d_sat->addClause(~d_encoding->perm(k, j));

    // handle the case of unique original row
    if (info.used == 1) {
        if (d_fixed.set(rows.back(), current_row)) {
            comment(2) << "[inv] " << rows.back() << " fixed to " << current_row
                       << endl;
            d_statistics.uniqueInv->inc();
            if (d_colors)
                allowonly(rows.back(), current_row);
        }
    }
}

bool LexminSolver::test_sat(const std::pair<size_t, size_t> &cell,
                            const std::vector<size_t> &vals) {
    const auto start_time = read_cpu_time();

    Minisat::vec<Minisat::Lit> assumps(1, mkLit(d_sat->fresh()));
    const auto &selector = assumps[0];
    d_encoding->encode_shot(cell, vals, selector);
    const auto rv = d_sat->solve(assumps);

    d_statistics.satTime->inc(read_cpu_time() - start_time);
    d_statistics.satCalls->inc();
    if (d_options.last_solution && rv)
        make_last_permutation();
    /* d_sat->addClause(rv ? selector : ~selector); */
    /* d_sat->addClause(~selector); */
    d_sat->releaseVar(~selector);
    /* TRACE(ccomment(3) << " " << std::get<2>(d_assignments.back()) << ":"
     */
    /*                   << SHOW_TIME(read_cpu_time() - start_time);); */
    return rv;
}
bool LexminSolver::test_sat(const Encoding::Assignment &asg) {
    const auto start_time = read_cpu_time();
    const auto rv = test_sat_inc(asg);
    d_statistics.satTime->inc(read_cpu_time() - start_time);
    d_statistics.satCalls->inc();
    if (d_options.last_solution && rv)
        make_last_permutation();
    TRACE(ccomment(3) << " " << std::get<2>(asg) << ":"
                      << SHOW_TIME(read_cpu_time() - start_time););
    return rv;
}

void LexminSolver::enc_inv_ord() {
    const auto n = d_table.order();
    BigInvariantCalculator big(d_output, d_table);
    big.calculate();
    const auto &invariants = big.invariants();
    std::vector<InvariantVector> sorted_invariants;
    InvariantVectorCmp cmp;
    for (const auto &[inv, info] : invariants)
        sorted_invariants.push_back(inv);
    std::sort(sorted_invariants.begin(), sorted_invariants.end(), cmp);

    std::vector<std::set<size_t>> blocks;
    for (const auto &inv : sorted_invariants)
        blocks.push_back(invariants.at(inv).elems);

    size_t range_start = 0;
    for (const auto &block : blocks) {
        const auto range_stop = range_start + block.size();
        if (d_options.verbose > 2) {
            comment(3) << "inv block " << range_start << ".."
                       << (range_stop - 1) << " {";
            for (const auto k : block)
                d_output.ccomment(3) << " " << k;
            d_output.ccomment(3) << " }" << std::endl;
        }
        for (const auto elem : block) {
            for (size_t i = 0; i < range_start; i++)
                d_sat->addClause(~d_encoding->perm(elem, i));
            for (size_t i = range_stop; i < n; i++)
                d_sat->addClause(~d_encoding->perm(elem, i));
        }
        if (block.size() == 1) {
            const auto elem = *block.begin();
            if (d_fixed.set(elem, range_start)) {
                comment(2) << elem << " fixed to " << d_fixed.src2dst(elem)
                           << " (inv_ord)" << std::endl;
                if (d_colors)
                    allowonly(elem, range_start);
            }
        }
        range_start = range_stop;
    }
    closure_fixed();
}

void LexminSolver::opt1stRow() {
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

    for (size_t i = 0; i < max_repeats; i++) {
        d_fixed_cells->set(0, i, 0);
        d_encoding->encode_pos({0, i, 0}, SATSPC::lit_Undef);
    }

    size_t some_first_row = -1; // some row can be 1st
    size_t count_can_be_first = 0;
    for (auto row = n; row--;) {
        const bool can_be_first = repeats[row] == max_repeats;
        if (can_be_first) {
            count_can_be_first++;
            some_first_row = row;
        } else {
            d_sat->addClause(~d_encoding->perm(row, 0));
            if (d_colors)
                try_disallow(row, 0);
        }
    }

    assert(count_can_be_first > 0);
    comment(3) << "count_can_be_first " << count_can_be_first << std::endl;

    // handle the case of unique first row
    if (count_can_be_first == 1) {
        d_statistics.unique1stRow->inc();
        d_fixed.set(some_first_row, 0);
        comment(2) << some_first_row << " fixed to 0 (opt1stRow)" << std::endl;
        if (d_colors)
            allowonly(some_first_row, 0);
        if (some_first_row == 0)
            d_statistics.unique1stRowFirst->inc();
    }
}

bool LexminSolver::test_sat_inc(const Encoding::Assignment &asg) {
    Minisat::vec<Minisat::Lit> assumps(1, mkLit(d_sat->fresh()));
    const auto &selector = assumps[0];
    d_encoding->encode_pos(asg, selector);
    const auto res = d_sat->solve(assumps);
    d_sat->addClause(res ? selector : ~selector);
    return res;
}

void LexminSolver::calculate_budgets_col() {
    assert(d_budgets);
    const auto n = d_table.order();
    size_t max_col_occs = 0; // max occurrences for non-fixed elements
    std::vector<size_t> max_col_occs_fixed(n, 0); // max occs. fixed elems
    std::vector<size_t> col_occs(n);              // occs. in current col
    for (auto col = n; col--;) {
        col_occs.assign(n, 0);
        for (auto row = n; row--;)
            col_occs[d_table.get(row, col)]++;

        // update occurrence maxima
        for (auto val = n; val--;) {
            auto &mx = d_fixed.is_fixed_src(val) ? max_col_occs_fixed[val]
                                                 : max_col_occs;
            mx = std::max(mx, col_occs[val]);
        }
    }

    comment(2) << "max occ./col: " << max_col_occs << std::endl;

    // set up budgets
    d_budgets->d_col_budget.assign(n, max_col_occs); // default
    for (auto original_val = n; original_val--;) {
        if (!d_fixed.is_fixed_src(original_val))
            continue;
        const auto image = d_fixed.src2dst(original_val);
        d_budgets->d_col_budget[image] = max_col_occs_fixed[original_val];
        comment(2) << "max occ./col for " << image << ": "
                   << d_budgets->d_col_budget[image] << std::endl;
    }
}

struct BudgetAux {
    BudgetAux(size_t n) : max_row_occs(0), max_row_occs_fixed(n, 0) {}
    size_t max_row_occs; // max occurrences for non-fixed elements
    std::vector<size_t> max_row_occs_fixed; // max occs. fixed elems
    void update_max(bool fixed, size_t val, size_t occs) {
        auto &mx = fixed ? max_row_occs_fixed[val] : max_row_occs;
        mx = std::max(mx, occs);
    }
};

void LexminSolver::calculate_budgets_row_tot() {
    const auto n = d_table.order();
    std::vector<size_t> total_occs(n, 0); // occs. in table per element
    std::vector<size_t> cur_row_occs;     // occs. in current row
    BudgetAux all(n), idems(n), noidems(n);
    for (auto row = n; row--;) {
        if (d_used[row]) // skip used rows
            continue;

        // update occurrences for each value in the row
        cur_row_occs.assign(n, 0);
        for (auto col = n; col--;) {
            const auto val = d_table.get(row, col);
            cur_row_occs[val]++;
            total_occs[val]++;
        }
        const bool use_idem =
            d_options.budget_idem && d_table.get(row, row) == row;
        const bool use_nonidem = d_options.budget_idem && !use_idem;

        // update row occurrence maxima
        for (auto val = n; val--;) {
            const auto occs = cur_row_occs[val];
            const auto fixed = d_fixed.is_fixed_src(val);
            all.update_max(fixed, val, occs);
            if (use_idem)
                idems.update_max(fixed, val, occs);
            if (use_nonidem)
                noidems.update_max(fixed, val, occs);
        }
    }

    // maximum number of total occurrences over non-fixed all elements
    size_t max_total_occs = 0;
    for (auto val = n; val--;)
        if (!d_fixed.is_fixed_src(val))
            max_total_occs = std::max(max_total_occs, total_occs[val]);

    TRACE(print_set(comment(4) << "total_occs ", total_occs) << std::endl;);
    comment(2) << "max all occ./row: " << all.max_row_occs << std::endl;
    comment(2) << "max idem occ./row: " << idems.max_row_occs << std::endl;
    comment(2) << "max noidem occ./row: " << noidems.max_row_occs << std::endl;
    comment(2) << "max occ./tot: " << max_total_occs << std::endl;

    // set up budgets
    d_budgets->d_all_budget.assign(n, all.max_row_occs);        // default
    d_budgets->d_idem_budget.assign(n, idems.max_row_occs);     // default
    d_budgets->d_noidem_budget.assign(n, noidems.max_row_occs); // default
    d_budgets->d_total_budget.assign(n, max_total_occs);        // default
    // set up budgets for fixed elements
    for (auto original_val = n; original_val--;) {
        if (!d_fixed.is_fixed_src(original_val))
            continue;
        const auto image = d_fixed.src2dst(original_val);
        d_budgets->d_all_budget[image] = all.max_row_occs_fixed[original_val];
        d_budgets->d_idem_budget[image] =
            idems.max_row_occs_fixed[original_val];
        d_budgets->d_noidem_budget[image] =
            noidems.max_row_occs_fixed[original_val];
        d_budgets->d_total_budget[image] = total_occs[original_val];
        /* comment(2) << "max occ./row for " << image << ": " */
        /*            << d_budgets->d_row_budget[image] << std::endl; */
        /* comment(2) << "max occ./tot for " << image << ": " */
        /*            << d_budgets->d_total_budget[image] << std::endl; */
    }
}

BinaryFunction *LexminSolver::make_solution() {
    assert(d_is_solved);
    const auto n = d_table.order();
    auto solution = new BinaryFunction(n);
    for (const auto &[row, col, val] : d_assignments)
        solution->set(row, col, val);
    return solution;
}

CompFunction LexminSolver::make_solution_comp() {
    const auto n = d_table.order();
    CompFunctionBuilder b(n, 2);
    [[maybe_unused]] size_t pos = 0;
    for (const auto &[row, col, val] : d_assignments) {
        assert(row == pos / n && col == pos % n);
        b.push(val);
        pos++;
    }
    return b.make();
}
