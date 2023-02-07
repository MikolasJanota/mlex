/*
 * File:  invariants.cpp
 * Author:  mikolas
 * Created on:  Mon Jan 9 15:50:15 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "invariants.h"
#include "binary_function.h"
#include "immutable_vector.h"
#include "options.h"
#include <cassert>
#include <cstddef>
#include <iostream> // for operator<<, ostream, basic_ostream, endl
#include <limits>
#include <vector>

/* #define INVARIANT_TRACING */
#if !defined(NDEBUG) || defined(INVARIANT_TRACING)
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

void Invariants::calculate() {
    const auto n = d_table.order();
    InvariantCalculator calc(n);

    for (auto row = n; row--;) {
        calc.set_row(row);
        for (auto col = n; col--;)
            calc.set_val(col, d_table.get(row, col));
        InvariantVector iinv = calc.make_ivec();
        auto [it, _] = d_invariants.insert({iinv, Info()});
        it->second.original_rows.push_back(row);
    }
    if (d_output.d_options.verbose > 2) {
        for (const auto &[inv, info] : d_invariants) {
            d_output.comment(3) << "inv " << inv << " {";
            for (const auto k : info.original_rows)
                d_output.ccomment(3) << " " << k;
            d_output.ccomment(3) << " }" << std::endl;
        }
    }
}
class Looping {
  public:
    Looping(Output &output, const std::vector<size_t> &fun)
        : d_output(output), d_order(fun.size()), d_fun(fun),
          d_value(d_order, std::numeric_limits<std::size_t>::max()){};

    bool has_val(size_t i) { return d_value[i] <= d_order; }

    size_t calc_loop(size_t query_ix) {
        if (has_val(query_ix)) {
            TRACE(d_output.comment(4)
                      << "lc: " << query_ix << ":" << d_value[query_ix]
                      << " (mem)" << std::endl;);
            return d_value[query_ix];
        }

        std::vector<size_t> time(d_order,
                                 std::numeric_limits<std::size_t>::max());
        std::vector<size_t> stack;
        auto next = query_ix;
        size_t t = 0;
        // either hit a cycle or something known
        while (!has_val(next) && time[next] >= d_order) {
            stack.push_back(next);
            time[next] = t++;
            next = d_fun[next];
        }

        const size_t known_sz =
            has_val(next) ? d_value[next] : (t - time[next]);

        if (!has_val(next)) { // process discovered cycle
            assert(time[next] < d_order);
            for (auto _ = known_sz; _--;) {
                d_value[stack.back()] = known_sz;
                stack.pop_back();
            }
        }

        // tail of the known part
        size_t edges = 1;
        while (!stack.empty()) {
            d_value[stack.back()] = known_sz + edges++;
            stack.pop_back();
        }

        TRACE(d_output.comment(4) << "lc: " << query_ix << ":"
                                  << d_value[query_ix] << std::endl;);
        assert(has_val(query_ix));
        return d_value[query_ix];
    }

  private:
    Output &d_output;
    const size_t d_order;
    const std::vector<size_t> &d_fun;
    std::vector<size_t> d_value;
};

void DiagInvariants::calculate() {
    using vec = std::vector<size_t>;
    std::vector<vec> invs(d_order, vec{0, 0});
    Looping lc(d_output, d_diagonal);
    for (auto i = d_order; i--;) {
        invs[d_diagonal[i]]
            [InvariantType::REPEATS]++; // number of occs of element
        invs[i][InvariantType::LOOP] = lc.calc_loop(i);
    }
    d_invariants.resize(d_order);
    for (auto i = d_order; i--;) {
        d_invariants[i] = InvariantVector(invs[i]);
        TRACE(d_output.comment(3)
                  << "diag inv: " << i << ":" << d_invariants[i] << std::endl;);
    }
}

void DiagInvariants::calc_inverse() {
    d_inv2elems = std::make_unique<inv_map>();
    for (auto i = d_order; i--;) {
        auto [it, _] = d_inv2elems->insert({d_invariants[i], Info()});
        it->second.elems.insert(i);
    }
}

void DiagInvariants::set(size_t i, size_t val) { d_diagonal[i] = val; }
