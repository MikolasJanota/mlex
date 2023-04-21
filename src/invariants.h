/*
 * File:  invariants.h
 * Author:  mikolas
 * Created on:  Mon Jan 9 15:38:58 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "auxiliary.h"
#include "immutable_vector.h"
#include "options.h"
#include <cassert>
#include <cstddef>
#include <list>
#include <limits>
#include <memory> // for unique_ptr
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility> // for pair
#include <vector>

typedef ImmutableVector<size_t> InvariantVector;

class BinaryFunction;
class Output;

class InvariantCalculator {
  public:
    InvariantCalculator(size_t n) : d_n(n) {}

    void add_loop(size_t loop_sz) {
        const auto inv_id = invariant_count + loop_sz;
        if (d_invv.size() <= inv_id)
            d_invv.resize(inv_id + 1, 0);
        d_invv[inv_id]++;
    }

    void set_val(size_t col, size_t val) {
        size_t inv_id = 0;
        if (val == d_row)
            d_invv[inv_id]++;
        inv_id++;
        if (val == col)
            d_invv[inv_id]++;
        inv_id++;
        if ((val == col) && (val == d_row))
            d_invv[inv_id]++;
        inv_id++;
        if (!d_seen[val]) {
            d_seen[val] = true;
            d_invv[inv_id]++;
        }
        inv_id++;
        assert(inv_id == invariant_count);
    }

    inline InvariantVector make_ivec() { return InvariantVector(d_invv); }

    void set_row(size_t row) {
        d_row = row;
        d_seen.clear();
        d_seen.resize(d_n, false);
        d_invv.clear();
        d_invv.resize(invariant_count, 0);
    }

    const size_t invariant_count = 4;

  private:
    const size_t d_n;
    size_t d_row;
    std::vector<size_t> d_invv;
    std::vector<bool> d_seen;
};

class DiagInvariants {
  public:
    enum InvariantType { REPEATS = 0, LOOP = 1 };
    struct Info {
        std::set<size_t> elems;
    };

    DiagInvariants(Output &output, size_t order)
        : d_output(output), d_order(order), d_diagonal(d_order, -1) {}

    void calculate();
    void calc_inverse();
    void set(size_t i, size_t val);

    void add(size_t i) {
        assert(i < d_order);
        d_elems.insert(i);
    }

    InvariantVector get_invariant(size_t i) const {
        assert(contains(d_elems, i));
        return d_invariants[i];
    }

    size_t get_reps(size_t i) const {
        assert(contains(d_elems, i));
        return d_invariants[i][InvariantType::REPEATS];
    }

    size_t get_loop(size_t i) const {
        assert(contains(d_elems, i));
        return d_invariants[i][InvariantType::LOOP];
    }

    const std::set<size_t> &get_elems(const InvariantVector &inv) const {
        return d_inv2elems->at(inv).elems;
    }

    const std::unordered_set<size_t> &elems() const { return d_elems; }

  private:
    Output &d_output;
    const size_t d_order;
    std::vector<size_t> d_diagonal;
    std::unordered_set<size_t> d_elems;
    std::vector<InvariantVector> d_invariants;
    using inv_map =
        std::unordered_map<InvariantVector, Info, ImmutableVector_hash<size_t>,
                           ImmutableVector_equal<size_t>>;
    std::unique_ptr<inv_map> d_inv2elems;
};

class Invariants {
  public:
    struct Info {
        std::list<size_t> original_rows;
        size_t used = 0;
    };

    Invariants(Output &output, const BinaryFunction &table)
        : d_output(output), d_table(table) {}

    void calculate();

    Info &get(const InvariantVector &invv) {
        const auto i = d_invariants.find(invv);
        assert(d_invariants.end() != i);
        return i->second;
    }

  private:
    Output &d_output;
    const BinaryFunction &d_table;
    std::unordered_map<InvariantVector, Info, ImmutableVector_hash<size_t>,
                       ImmutableVector_equal<size_t>>
        d_invariants;
};

class Looping {
  public:
    Looping(Output &output, const std::vector<size_t> &fun)
        : d_output(output), d_order(fun.size()), d_fun(fun),
          d_value(d_order, std::numeric_limits<std::size_t>::max()){};

    bool has_val(size_t i) { return d_value[i] <= d_order; }
    size_t calc_loop(size_t query_ix);

  private:
    Output &d_output;
    const size_t d_order;
    const std::vector<size_t> &d_fun;
    std::vector<size_t> d_value;
};

