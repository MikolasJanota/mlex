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
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
using SATSPC::Lit;
using SATSPC::mkLit;

#define VERB(level, code)                                                      \
    if (d_options.verbose >= level)                                            \
        do {                                                                   \
            code                                                               \
    } while (0)

#ifdef NDEBUG
#define TRACE(code)
#else
#define TRACE(code)                                                            \
    do {                                                                       \
        code                                                                   \
    } while (0);                                                               \
    do {                                                                       \
        std::cout.flush();                                                     \
    } while (0)
#endif

class Budgets {
  public:
    Budgets(const BinaryFunction &table) : d_table(table) {}

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

LexminSolver::LexminSolver(Output &output, const BinaryFunction &table)
    : d_output(output), d_options(output.d_options),
      d_statistics(output.d_statistics), d_table(table),
      d_invariants(output, table) {}
LexminSolver::~LexminSolver() {}

size_t LexminSolver::find_value(const std::optional<size_t> &last_val) {
    switch (d_options.search_type) {
    case SearchType::lin_us: return find_value_unsat_sat(last_val);
    case SearchType::lin_su: return find_value_sat_unsat(last_val);
    case SearchType::bin: return find_value_bin(last_val);
    }
    assert(false);
    return -1;
}

size_t LexminSolver::find_value_bin(const std::optional<size_t> &last_val) {
    const auto n = d_table.order();
    auto &[row, col, cur_val] = d_assignments.back();
    const auto cell = std::make_pair<>(row, col);

    TRACE(comment(3) << "(" << row << " " << col << ") :";);
    size_t ub = last_val ? (*last_val + 1) : n; // upper bound
    TRACE(comment(3) << "(iub:" << ub << ") ";);
    std::vector<size_t> a, b;
    std::vector<size_t> *vals = &a, *top = &b;

    for (size_t v = 0; v < ub; v++) {
        if (!d_budgets || d_budgets->has_budget(col, v))
            vals->push_back(v);
    }

    while (vals->size() > 1) {
        VERB(4, print_set(comment(4) << "vals ", *vals););

        assert(top->empty());
        const auto split = vals->size() / 2;
        for (size_t h = split; h < vals->size(); h++)
            top->push_back(vals->at(h));
        vals->resize(split);

        if (test_sat(cell, *vals)) {
            assert(!d_last_permutation.empty());
            ub = std::min(ub, get_val(row, col));
            TRACE(comment(3) << "(ub:" << ub << ") ";);
            while (!vals->empty() && vals->back() > ub)
                vals->pop_back();
        } else {
            std::swap(vals, top);
        }
        top->clear();
    }

    assert(vals->size() == 1);
    cur_val = vals->back();
    VERB(4, comment(4) << "val: " << cur_val;);
    d_encoding->encode_pos(d_assignments.back(), SATSPC::lit_Undef);
    return cur_val;
}

size_t
LexminSolver::find_value_sat_unsat(const std::optional<size_t> &last_val) {
    const auto n = d_table.order();
    auto &[row, col, cur_val] = d_assignments.back();
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
    VERB(4, comment(4) << "val: " << cur_val;);
    d_encoding->encode_pos(d_assignments.back(), SATSPC::lit_Undef);
    return ub;
}

size_t
LexminSolver::find_value_unsat_sat(const std::optional<size_t> &last_val) {
    auto &[row, col, cur_val] = d_assignments.back();

    TRACE(comment(3) << "(" << row << " " << col << ") :";);

    for (;; cur_val++) {
        assert(cur_val < d_table.order());
        TRACE(if (d_budgets) d_budgets->print(ccomment(4), col, cur_val)
                  << " ";);

        const bool skip = d_budgets && !d_budgets->has_budget(col, cur_val);
        auto found = !skip && last_val && *last_val == cur_val;
        if (found) {
            TRACE(ccomment(3) << " " << cur_val << ":FREE";);
            if (d_options.incremental)
                d_encoding->encode_pos(d_assignments.back(), SATSPC::lit_Undef);

        } else {
            found = !skip && test_sat();
        }

        if (d_budgets && found)
            d_budgets->dec_budget(col, cur_val);

        if (found)
            return cur_val;
    }
    assert(false);
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

    if (budgeting) {
        d_budgets = std::make_unique<Budgets>(d_table);
        calculate_budgets_row_tot();
        calculate_budgets_col();
        d_budgets->reset_cols_budget();
    }

    d_assignments.clear();
    d_assignments.reserve(n * n);

    InvariantCalculator calc(n);

    for (size_t row = 0; row < n; row++) {
        if (budgeting)
            d_budgets->reset_cur_row_budget();

        if (d_options.invariants)
            calc.set_row(row);

        for (size_t col = 0; col < n; col++) {
            d_assignments.push_back({row, col, 0});
            const auto cur_val =
                find_value(d_last_permutation.empty()
                               ? std::nullopt
                               : std::make_optional(get_val(row, col)));

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

    make_solution();
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

    VERB(3, print_set(comment(3) << "used ", rows););

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
    /* TRACE(ccomment(3) << " " << std::get<2>(d_assignments.back()) << ":" */
    /*                   << SHOW_TIME(read_cpu_time() - start_time);); */
    return rv;
}
bool LexminSolver::test_sat() {
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
    if (d_options.opt1stRow)
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

    VERB(3, print_set(comment(3) << "idems ", idems) << std::endl;);

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

    VERB(3, print_set(comment(3) << "total_occs ", total_occs) << std::endl;);
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

void LexminSolver::make_solution() {
    const auto n = d_table.order();
    d_solution = std::make_unique<BinaryFunction>(n);
    for (const auto &[row, col, val] : d_assignments)
        d_solution->set(row, col, val);
}

void LexminSolver::print_solution(std::ostream &output) {
    if (d_options.mace_format)
        d_solution->print_mace(output);
    else
        d_solution->print_gap(output);
}

