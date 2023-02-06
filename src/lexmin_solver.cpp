/*
 * File:  lexmin_solver.cpp
 * Author:  mikolas
 * Created on:  Wed Dec 14 13:29:54 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "lexmin_solver.h"
#include "auxiliary.h"
#include "encoding.h"
#include "invariants.h"
#include "minisat/core/SolverTypes.h"
#include "minisat_ext.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <tuple> // for get
#include <vector>
#include <zlib.h>
using SATSPC::Lit;
using SATSPC::mkLit;

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

    virtual bool has_budget(size_t val) const {
        return (val != d_row || d_idem_count) && d_budgets[val];
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

  private:
    size_t d_row;
    std::vector<size_t> d_budgets;
    size_t d_idem_count;
};

class RowBudgets {
  public:
    RowBudgets(const BinaryFunction &table) : d_table(table) {}

    bool has_budget(size_t col, size_t val) const {
        assert(col < d_table.order());
        assert(val < d_table.order());
        return (d_current_row_budget[val] > 0 && d_total_budget[val] > 0 &&
                d_cols_budget[col][val] > 0);
    }

    void dec_budget(size_t col, size_t val) {
        assert(col < d_table.order());
        assert(val < d_table.order());
        d_current_row_budget[val]--;
        d_total_budget[val]--;
        d_cols_budget[col][val]--;
    }

    std::ostream &print(std::ostream &out, size_t col, size_t val) const {
        return out << "(budg r/c/t: " << d_current_row_budget[val] << " "
                   << d_cols_budget[col][val] << " " << d_total_budget[val]
                   << ")";
    }

    void reset_cur_row_budget() {
        assert(d_row_budget.size() == d_table.order());
        d_current_row_budget = d_row_budget;
    }

    void reset_cols_budget() {
        assert(d_col_budget.size() == d_table.order());
        d_cols_budget.assign(d_table.order(), d_col_budget);
    }

    std::vector<std::vector<size_t>> d_cols_budget;
    std::vector<size_t> d_row_budget;
    std::vector<size_t> d_col_budget;
    std::vector<size_t> d_total_budget;
    std::vector<size_t> d_current_row_budget;

  private:
    const BinaryFunction &d_table;
};

class RowBudgeter : public IBudget {
  public:
    RowBudgeter(RowBudgets *paren, bool trivial)
        : d_paren(paren), d_triv(trivial), d_col(0) {}

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
    : d_output(output), d_options(output.d_options),
      d_statistics(output.d_statistics), d_table(table),
      d_invariants(output, table) {}

LexminSolver::~LexminSolver() {}

size_t LexminSolver::find_value(Encoding::Assignment &asg, IBudget &budget,
                                const std::optional<size_t> &last_val) {
    if (d_fixed_values) {
        auto &[row, col, val] = asg;
        if (d_fixed_values->is_set(row, col)) {
            val = d_fixed_values->get(row, col);
            TRACE(comment(3)
                      << "fixed cell: (" << row << "," << col << ")=" << val;);
            return val;
        }
    }

    switch (d_options.search_type) {
    case SearchType::lin_us: return find_value_unsat_sat(asg, last_val);
    case SearchType::lin_su: return find_value_sat_unsat(asg, last_val);
    case SearchType::bin: return find_value_bin(asg, last_val);
    case SearchType::bin2: return find_value_bin2(asg, budget, last_val);
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
        if (!d_budgets || budget.has_budget(v))
            vals->push_back(v);
    }

    while (!vals->empty()) {
        assert(top->empty());
        TRACE(print_set(comment(4) << "vals ", *vals););

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

size_t LexminSolver::find_value_bin(Encoding::Assignment &asg,
                                    const std::optional<size_t> &last_val) {
    assert(false); // TODO
    const auto n = d_table.order();
    auto &[row, col, cur_val] = asg;
    const auto cell = std::make_pair<>(row, col);

    TRACE(comment(3) << "(" << row << " " << col << ") :";);
    size_t ub = last_val ? *last_val : n; // upper bound
    TRACE(comment(3) << "(iub:" << ub << ") ";);
    std::vector<size_t> a, b;
    std::vector<size_t> *vals = &a, *top = &b;

    for (size_t v = 0; v < ub; v++) {
        if (!d_budgets || d_budgets->has_budget(col, v))
            vals->push_back(v);
    }

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
LexminSolver::find_value_sat_unsat(Encoding::Assignment &asg,
                                   const std::optional<size_t> &last_val) {

    assert(false); // TODO
    const auto n = d_table.order();
    auto &[row, col, cur_val] = asg;
    const auto cell = std::make_pair<>(row, col);

    TRACE(comment(3) << "(" << row << " " << col << ") :";);
    size_t ub = last_val ? *last_val : n; // upper bound
    TRACE(comment(3) << "(iub:" << ub << ") ";);
    std::vector<size_t> vals;
    for (size_t v = 0; v < ub; v++) {
        if (!d_budgets || d_budgets->has_budget(col, v))
            vals.push_back(v);
    }
    /* print_set(std::cerr << "in:", vals) << std::endl; */
    while (!vals.empty() && test_sat(cell, vals)) {
        /* print_set(std::cerr << "bf:", vals) << std::endl; */
        assert(!d_last_permutation.empty());
        ub = std::min(ub, get_val(row, col));
        TRACE(comment(3) << "(ub:" << ub << ") ";);
        while (!vals.empty() && vals.back() >= ub)
            vals.pop_back();
        /* print_set(std::cerr << "af:", vals) << std::endl; */
    }
    cur_val = ub;
    TRACE(comment(4) << "val: " << cur_val;);
    d_encoding->encode_pos(asg, SATSPC::lit_Undef);
    return ub;
}

size_t
LexminSolver::find_value_unsat_sat(Encoding::Assignment &asg,
                                   const std::optional<size_t> &last_val) {
    assert(false); // TODO
    auto &[row, col, cur_val] = asg;

    TRACE(comment(3) << "(" << row << " " << col << ") :";);

    for (;; cur_val++) {
        assert(cur_val < d_table.order());
        TRACE(if (d_budgets) d_budgets->print(ccomment(5), col, cur_val)
                  << " ";);

        const bool skip = d_budgets && !d_budgets->has_budget(col, cur_val);
        auto found = !skip && last_val && *last_val == cur_val;
        if (found) {
            TRACE(ccomment(3) << " " << cur_val << ":FREE";);
            if (d_options.incremental)
                d_encoding->encode_pos(asg, SATSPC::lit_Undef);

        } else {
            found = !skip && test_sat();
        }

        if (found)
            return cur_val;
    }
    assert(false);
}

void LexminSolver::run_diagonal() {
    using std::make_pair;
    using std::max;
    using std::optional;
    using std::pair;
    const auto n = d_table.order();

    if (!d_fixed_values)
        d_fixed_values = std::make_unique<BinaryFunction>(n);

    DiagInvariants di_orig(d_output, n);
    for (size_t i = 0; i < n; i++)
        di_orig.set(i, d_table.get(i, i));
    di_orig.calculate();
    size_t max_idem_reps = 0;
    size_t max_reps = 0;
    std::vector<size_t> idems;
    for (size_t i = 0; i < n; i++) {
        const auto &inv = di_orig.get_invariant(i);
        const auto reps = inv[DiagInvariants::InvariantType::REPEATS];
        max_reps = max(max_reps, reps);
        if (inv[DiagInvariants::InvariantType::LOOP] == 1) {
            // idempotent
            max_idem_reps = max(max_idem_reps, reps);
            idems.push_back(i);
        }
    }
    if (!idems.empty()) {
        std::set<size_t> max_rep_idems;
        for (auto i : idems)
            if (di_orig.get_invariant(
                    i)[DiagInvariants::InvariantType::REPEATS] == max_idem_reps)
                max_rep_idems.insert(i);
        if (max_rep_idems.size() == 1) {
            d_statistics.uniqueDiag1Elem->inc();
            const auto i = *(max_rep_idems.begin());
            d_fixed[i] = 0;
            comment(2) << i << " fixed to 0 (diagonal)" << std::endl;
        }
        for (auto i = n; i--;)
            if (!contains(max_rep_idems, i))
                d_sat->addClause(~d_encoding->perm(i, 0));
    }

    /* comment(2) << "max idemp diag reps: " << max_idem_reps << std::endl; */
    /* comment(2) << "max noidemp diag reps: " << max_no_idem_reps << std::endl;
     */
    DiagBudget budg(n, max_reps, idems.size());
    for (size_t i = 0; i < n; i++) {
        Encoding::Assignment asg{i, i, 0};
        const auto last_val = d_last_permutation.empty()
                                  ? std::nullopt
                                  : std::make_optional(get_val(i, i));
        budg.set_row(i);
        const auto new_val = find_value(asg, budg, last_val);
        d_fixed_values->set(i, i, new_val);
        budg.dec(new_val);
        TRACE(ccomment(3) << std::endl;);
    }
    d_statistics.satCalls->print(comment(2) << "diag:") << std::endl;

    DiagInvariants di_new(d_output, n);
    for (size_t i = 0; i < n; i++)
        di_new.set(i, d_fixed_values->get(i, i));
    di_new.calculate();
    di_new.calc_inverse();

    std::set<size_t> fixed_elements;
    for (size_t i = 0; i < n; i++) {
        const auto inv = di_orig.get_invariant(i);
        const auto &inv_pre_image = di_new.get_info(inv).original_elems;
        // disable permuting to nonmatching elements
        for (size_t j = 0; j < n; j++)
            if (!contains(inv_pre_image, j))
                d_sat->addClause(~d_encoding->perm(i, j));
        // handle the case of unique pre-image
        if (inv_pre_image.size() == 1) {
            d_statistics.uniqueDiagElem->inc();
            const auto uniq_new = *(inv_pre_image.begin());
            d_fixed[i] = uniq_new;
            fixed_elements.insert(i);
            comment(2) << i << " fixed to " << uniq_new << " (diag)"
                       << std::endl;
        }
    }
    for (auto row : fixed_elements)
        for (auto col : fixed_elements) {
            const auto orig_val = d_table.get(row, col);
            if (contains(fixed_elements, orig_val)) {
                d_fixed_values->set(d_fixed[row], d_fixed[col],
                                    d_fixed[orig_val]);
                comment(2) << "fixed cell " << row << " " << col << std::endl;
            }
        }
}

void LexminSolver::solve() {
    const auto n = d_table.order();
    const auto budgeting = d_options.budgeting;
    d_fixed.resize(n, n);
    d_used.resize(n, false);

    if (d_options.invariants) {
        d_invariants.calculate();
    }

    if (d_options.incremental) {
        d_sat = std::make_unique<SATSPC::MiniSatExt>();
        make_encoding();
        d_encoding->encode_bij();
    }

    if (d_options.diagonal) {
        run_diagonal();
    }

    if (budgeting) {
        d_budgets = std::make_unique<RowBudgets>(d_table);
        calculate_budgets_row_tot();
        calculate_budgets_col();
        d_budgets->reset_cols_budget();
    }

    d_assignments.clear();
    d_assignments.reserve(n * n);

    InvariantCalculator calc(n);

    const auto start_time = read_cpu_time();
    for (size_t row = 0; row < n; row++) {
        comment(3) << "row:" << row << " "
                   << SHOW_TIME(read_cpu_time() - start_time) << std::endl;

        if (budgeting)
            d_budgets->reset_cur_row_budget();

        if (d_options.invariants)
            calc.set_row(row);

        RowBudgeter budgeter(d_budgets.get(), !d_budgets);
        for (size_t col = 0; col < n; col++) {
            d_assignments.push_back({row, col, 0});
            find_value(d_assignments.back(), budgeter,
                       d_last_permutation.empty()
                           ? std::nullopt
                           : std::make_optional(get_val(row, col)));
            const auto cur_val = std::get<2>(d_assignments.back());

            if (d_budgets)
                d_budgets->dec_budget(col, cur_val);

            if (d_options.invariants)
                calc.set_val(col, cur_val);

            TRACE(ccomment(3) << std::endl;);
        }

        bool update_budgets = false;

        if (d_options.invariants && row + 1 < n)
            update_budgets = process_invariant(calc.make_ivec(), row);

        if (row == 1 && d_0preimage) {
            d_used[*d_0preimage] = true; // mark preimage of first row as used
            update_budgets = true;
        }

        if (update_budgets && budgeting)
            calculate_budgets_row_tot();
    }
    d_is_solved = true;
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

    TRACE(print_set(comment(3) << "used ", rows););

    // disable further rows to be mapped to the used ones
    for (auto k : rows)
        for (size_t j = current_row + 1; j < d_table.order(); j++)
            d_sat->addClause(~d_encoding->perm(k, j));

    // handled the case of unique original row
    if (info.used == 1) {
        d_fixed[rows.back()] = current_row;
        comment(2) << rows.back() << " fixed to " << current_row << std::endl;
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
    d_sat->addClause(~selector);
    /* TRACE(ccomment(3) << " " << std::get<2>(d_assignments.back()) << ":"
     */
    /*                   << SHOW_TIME(read_cpu_time() - start_time);); */
    return rv;
}
bool LexminSolver::test_sat() {
    assert(false); // TODO
    const auto start_time = read_cpu_time();
    const auto rv = d_options.incremental ? test_sat_inc() : test_sat_noinc();
    d_statistics.satTime->inc(read_cpu_time() - start_time);
    d_statistics.satCalls->inc();
    if (d_options.last_solution && rv)
        make_last_permutation();
    TRACE(ccomment(3) << " " << std::get<2>(d_assignments.back()) << ":"
                      << SHOW_TIME(read_cpu_time() - start_time););
    return rv;
}

void LexminSolver::make_encoding() {
    d_encoding = std::make_unique<Encoding>(d_output, *d_sat, d_table);
    if (d_options.opt1stRow && !d_options.diagonal)
        opt1stRow();
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

    // identify a rows that can be first
    std::vector<bool> can_be_first(n, false);
    for (const auto row : idems)
        can_be_first[row] = repeats[row] == max_repeats;

    // process rows that can be first
    size_t maxRow = -1;
    size_t count_can_be_first = 0;
    for (auto row = n; row--;) {
        if (can_be_first[row]) {
            count_can_be_first++;
            maxRow = row;
        } else {
            d_sat->addClause(~d_encoding->perm(row, 0));
        }
    }

    assert(count_can_be_first);
    comment(3) << "count_can_be_first " << count_can_be_first << std::endl;

    // handle the case of unique first row
    if (count_can_be_first == 1) {
        d_statistics.unique1stRow->inc();
        d_fixed[maxRow] = 0;
        d_0preimage = maxRow;
        comment(2) << *d_0preimage << " fixed to 0 (opt1stRow)" << std::endl;
    }
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
            auto &mx = is_fixed(val) ? max_col_occs_fixed[val] : max_col_occs;
            mx = std::max(mx, col_occs[val]);
        }
    }

    comment(2) << "max occ./col: " << max_col_occs << std::endl;

    // set up budgets
    d_budgets->d_col_budget.assign(n, max_col_occs); // default
    for (auto original_val = n; original_val--;) {
        if (!is_fixed(original_val))
            continue;
        const auto image = d_fixed[original_val];
        d_budgets->d_col_budget[image] = max_col_occs_fixed[original_val];
        comment(2) << "max occ./col for " << image << ": "
                   << d_budgets->d_col_budget[image] << std::endl;
    }
}

void LexminSolver::calculate_budgets_row_tot() {
    const auto n = d_table.order();
    size_t max_row_occs = 0; // max occurrences for non-fixed elements
    std::vector<size_t> max_row_occs_fixed(n, 0); // max occs. fixed elems
    std::vector<size_t> row_occs(n, 0);           // occs. in current row
    std::vector<size_t> total_occs(n, 0);         // occs. in table per element
    for (auto row = n; row--;) {
        if (d_used[row]) // skip used rows
            continue;

        // update occurrences for each value in the row
        row_occs.assign(n, 0);
        for (auto col = n; col--;) {
            const auto val = d_table.get(row, col);
            row_occs[val]++;
            total_occs[val]++;
        }

        // update row occurrence maxima
        for (auto val = n; val--;) {
            auto &mx = is_fixed(val) ? max_row_occs_fixed[val] : max_row_occs;
            mx = std::max(mx, row_occs[val]);
        }
    }

    // maximum number of total occurrences over non-fixed all elements
    size_t max_total_occs = 0;
    for (auto val = n; val--;)
        if (!is_fixed(val))
            max_total_occs = std::max(max_total_occs, total_occs[val]);

    TRACE(print_set(comment(3) << "total_occs ", total_occs) << std::endl;);
    comment(2) << "max occ./row: " << max_row_occs << std::endl;
    comment(2) << "max occ./tot: " << max_total_occs << std::endl;

    // set up budgets
    d_budgets->d_row_budget.assign(n, max_row_occs);     // default
    d_budgets->d_total_budget.assign(n, max_total_occs); // default
    // set up budgets for fixed elements
    for (auto original_val = n; original_val--;) {
        if (!is_fixed(original_val))
            continue;
        const auto image = d_fixed[original_val];
        d_budgets->d_row_budget[image] = max_row_occs_fixed[original_val];
        d_budgets->d_total_budget[image] = total_occs[original_val];
        comment(2) << "max occ./row for " << image << ": "
                   << d_budgets->d_row_budget[image] << std::endl;
        comment(2) << "max occ./tot for " << image << ": "
                   << d_budgets->d_total_budget[image] << std::endl;
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
