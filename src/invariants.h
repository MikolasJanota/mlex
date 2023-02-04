/*
 * File:  invariants.h
 * Author:  mikolas
 * Created on:  Mon Jan 9 15:38:58 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "auxiliary.h"
#include "binary_function.h"
#include "immutable_vector.h"
#include "options.h"
#include <cstddef>
#include <list>
#include <set>
#include <unordered_map>
#include <vector>

typedef ImmutableVector<size_t> InvariantVector;

class InvariantCalculator {
  public:
    InvariantCalculator(size_t n) : d_n(n) {}

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
    struct Info {
        std::set<size_t> original_elems;
    };

    DiagInvariants(Output &output, size_t order)
        : d_output(output), d_options(output.d_options), d_order(order),
          d_diagonal(d_order, -1) {}

    void calculate();
    void calc_inverse();
    void set(size_t i, size_t val);
    InvariantVector get_invariant(size_t i) const { return d_invariants[i]; }
    Info get_info(const InvariantVector &inv) const {
        return d_inv2elems->at(inv);
    }

  private:
    Output &d_output;
    const Options &d_options;
    const size_t d_order;
    std::vector<size_t> d_diagonal;
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
        : d_output(output), d_options(output.d_options), d_table(table) {}

    void calculate();

    Info &get(const InvariantVector &invv) {
        const auto i = d_invariants.find(invv);
        assert(d_invariants.end() != i);
        return i->second;
    }

  private:
    Output &d_output;
    const Options &d_options;
    const BinaryFunction &d_table;
    std::unordered_map<InvariantVector, Info, ImmutableVector_hash<size_t>,
                       ImmutableVector_equal<size_t>>
        d_invariants;
};
