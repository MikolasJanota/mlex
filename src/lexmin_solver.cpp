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

    std::vector<std::vector<size_t>> cols_budget(n, d_col_budget);

    if (budgeting) {
        calculate_budgetsRowTot();
        calculate_budgetsCol();
        assert(d_col_budget.size() == n);
        cols_budget.assign(n, d_col_budget);
    }

    d_assignments.clear();
    d_assignments.reserve(n * n);

    InvariantCalculator calc(n);
    std::vector<size_t> current_row_budget;

    for (size_t row = 0; row < n; row++) {
        if (budgeting)
            current_row_budget = d_row_budget;

        if (d_options.invariants)
            calc.set_row(row);

        assert(!budgeting || d_total_budget.size() == n);
        assert(!budgeting || current_row_budget.size() == n);

        for (size_t col = 0; col < n; col++) {
            const bool has_last_eval = !d_last_permutation.empty();
            const auto last_val = has_last_eval ? get_val(row, col) : -1;
            assert(!budgeting || cols_budget[col].size() == n);
            d_assignments.push_back({row, col, 0});
            auto &cur_val = std::get<2>(d_assignments.back());

            TRACE(comment(3) << "(" << row << " " << col << ") :";);

            for (bool found = false; !found;) {
                TRACE(if (budgeting) ccomment(4)
                          << "(budg r/c/t: " << current_row_budget[cur_val]
                          << " " << cols_budget[col][cur_val] << " "
                          << d_total_budget[cur_val] << ") ";);

                const bool skip =
                    budgeting && (current_row_budget[cur_val] == 0 ||
                                  d_total_budget[cur_val] == 0 ||
                                  cols_budget[col][cur_val] == 0);
                found = !skip && has_last_eval && last_val == cur_val;
                if (found) {
                    TRACE(ccomment(3) << " " << cur_val << ":FREE";);
                    if (d_options.incremental)
                        d_encoding->encode_pos(d_assignments.back(),
                                               SATSPC::lit_Undef);

                } else {
                    found = !skip && test_sat();
                }

                if (budgeting && found) {
                    current_row_budget[cur_val]--;
                    d_total_budget[cur_val]--;
                    cols_budget[col][cur_val]--;
                }
                if (!found)
                    cur_val++;
            }
            assert(cur_val < n);

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

        if (update_budgets)
            calculate_budgetsRowTot();
    }
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

void LexminSolver::calculate_budgetsCol() {
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
    d_col_budget.assign(n, max_col_occs); // default
    for (auto original_val = n; original_val--;) {
        if (!is_fixed(original_val))
            continue;
        const auto image = d_fixed[original_val];
        d_col_budget[image] = max_col_occs_fixed[original_val];
        comment(2) << "max occ./col for " << image << ": "
                   << d_col_budget[image] << std::endl;
    }
}

void LexminSolver::calculate_budgetsRowTot() {
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

    // maximum number of total occurrences  over non-fixed all elements
    size_t max_total_occs = 0;
    for (auto val = n; val--;)
        if (!is_fixed(val))
            max_total_occs = std::max(max_total_occs, total_occs[val]);

    VERB(3, print_set(comment(3) << "total_occs ", total_occs) << std::endl;);
    comment(2) << "max occ./row: " << max_row_occs << std::endl;
    comment(2) << "max occ./tot: " << max_total_occs << std::endl;

    // set up budgets
    d_row_budget.assign(n, max_row_occs);     // default
    d_total_budget.assign(n, max_total_occs); // default
    for (auto original_val = n; original_val--;) {
        if (!is_fixed(original_val))
            continue;
        const auto image = d_fixed[original_val];
        d_row_budget[image] = max_row_occs_fixed[original_val];
        d_total_budget[image] = total_occs[original_val];
        comment(2) << "max occ./row for " << image << ": "
                   << d_row_budget[image] << std::endl;
        comment(2) << "max occ./tot for " << image << ": "
                   << d_total_budget[image] << std::endl;
    }
}

BinaryFunction *LexminSolver::make_solution() {
    const auto n = d_table.order();
    auto solution = new BinaryFunction(n);
    for (const auto &[row, col, val] : d_assignments)
        solution->set(row, col, val);
    return solution;
}

CompFunction LexminSolver::make_solution_comp() {
    const auto n = d_table.order();
    CompFunctionBuilder b(n, 2);
    size_t pos = 0;
    for (const auto &[row, col, val] : d_assignments) {
        assert(row == pos / n && col == pos % n);
        b.push(val);
        pos++;
    }
    return b.make(true);
}
