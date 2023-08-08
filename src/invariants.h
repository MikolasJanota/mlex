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
#include <cassert>
#include <cstddef>
#include <limits>
#include <list>
#include <memory> // for unique_ptr
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility> // for pair
#include <vector>

typedef ImmutableVector<size_t> InvariantVector;

class BinaryFunction;
class Output;

class BigInvariant {
  public:
    enum InvariantType {
        IDEM = 0,
        ROW_ORDER,
        ROW_DIST,
        COL_ORDER,
        COL_DIST,
        DIAG_FREQ,
        DIAG_ORDER,
        TOT_FREQ,
        LAST
    };

    BigInvariant(size_t order) : d_order(order) {
        const auto n = d_order;
        d_vals.resize(InvariantType::LAST);
        get_vals(InvariantType::IDEM).resize(1);
        get_vals(InvariantType::ROW_ORDER).resize(n + 1);
        get_vals(InvariantType::COL_ORDER).resize(n + 1);
        get_vals(InvariantType::ROW_DIST).resize(n + 1);
        get_vals(InvariantType::COL_DIST).resize(n + 1);
        get_vals(InvariantType::DIAG_FREQ).resize(1);
        get_vals(InvariantType::DIAG_ORDER).resize(1);
        get_vals(InvariantType::TOT_FREQ).resize(1);
    }

    std::vector<size_t> &get_vals(InvariantType it) {
        const auto ix = static_cast<size_t>(it);
        assert(ix < d_vals.size());
        return d_vals[ix];
    }

    inline InvariantVector make_ivec();
    std::ostream &print(std::ostream &o, const char *pref) {
        o << pref << "IDEM:" << get_vals(InvariantType::IDEM)[0] << std::endl;
        o << pref << "TOT_FREQ:" << get_vals(InvariantType::TOT_FREQ)[0]
          << std::endl;
        o << pref << "DIAG_ORDER:" << get_vals(InvariantType::DIAG_ORDER)[0]
          << std::endl;
        o << pref << "DIAG_FREQ:" << get_vals(InvariantType::DIAG_FREQ)[0]
          << std::endl;
        print_histo(o << pref << "ROW_ORDER:",
                    get_vals(InvariantType::ROW_ORDER))
            << std::endl;
        print_histo(o << pref << "COL_ORDER:",
                    get_vals(InvariantType::COL_ORDER))
            << std::endl;
        print_histo(o << pref << "ROW_DIST:", get_vals(InvariantType::ROW_DIST))
            << std::endl;
        print_histo(o << pref << "COL_DIST:", get_vals(InvariantType::COL_DIST))
            << std::endl;
        return o;
    }

  private:
    const size_t d_order;
    std::vector<std::vector<size_t>> d_vals;
    std::ostream &print_histo(std::ostream &o, const std::vector<size_t> &vs) {
        for (size_t i = 0; i < vs.size(); i++) {
            if (vs[i] > 0)
                o << " " << i << ":" << vs[i];
        }
        return o;
    }
};

/*  Utility class to produce invariants. */
class InvariantCalculator {
  public:
    InvariantCalculator(size_t n) : d_n(n) {}

    inline InvariantVector make_ivec() { return InvariantVector(d_invv); }

    /* resets the class for constructing invariant for a new row */
    void set_row(size_t row) {
        d_row = row;
        d_seen.clear();
        d_seen.resize(d_n, false);
        d_invv.clear();
        d_invv.resize(fixed_invariant_count, 0);
    }

    /* add looping value for an element (calculated by Looping) */
    void add_loop(size_t loop_sz) {
        assert(loop_sz <= d_n);
        const auto inv_id = fixed_invariant_count + loop_sz;
        if (d_invv.size() <= inv_id)
            d_invv.resize(inv_id + 1, 0);
        d_invv[inv_id]++;
    }

    /* add distance to row value for an element (calculated by Distances) */
    void add_distance(size_t distance) {
        assert(distance <= d_n);
        const auto inv_id = fixed_invariant_count + d_n + 1 + distance;
        if (d_invv.size() <= inv_id)
            d_invv.resize(inv_id + 1, 0);
        d_invv[inv_id]++;
    }

    /* register a value val in column col */
    void set_val(size_t col, size_t val) {
        size_t inv_id = 0;
        if ((val == col) && (val == d_row)) // idempotent
            d_invv[inv_id]++;
        inv_id++;
        if (val == d_row) // freq of r in row r
            d_invv[inv_id]++;
        inv_id++;
        // TODO: next useless with loops?
        if (val == col) // freq of t[r,c]=c
            d_invv[inv_id]++;
        inv_id++;
        if (!d_seen[val]) { // number of elements in the row
            d_seen[val] = true;
            d_invv[inv_id]++;
        }
        inv_id++;
        assert(inv_id == fixed_invariant_count);
    }

    static const size_t fixed_invariant_count = 4;

  private:
    const size_t d_n;
    size_t d_row;
    std::vector<size_t> d_invv;
    std::vector<bool> d_seen;
};

struct InvariantVectorCmp {
    inline bool operator()(const InvariantVector &v1,
                           const InvariantVector &v2) const {
        const auto fixedsz = InvariantCalculator::fixed_invariant_count;
        assert(v1.size() >= fixedsz);
        assert(v2.size() >= fixedsz);
        for (size_t i = 0; i < fixedsz; i++)
            if (v1[i] != v2[i])
                return v1[i] > v2[i];
        if (v1.size() != v2.size())
            return v1.size() > v2.size();
        for (size_t i = fixedsz; i < v1.size(); i++) {
            if (v1[i] != v2[i])
                return v1[i] > v2[i];
        }
        return false;
    }
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

/* The class maintains the set of invariants for the multiplication table.
 */
class Invariants {
  public:
    /* For invariant remember to which rows it belongs and how many times it
     * was already used in the target table. */
    struct Info {
        std::list<size_t> original_rows;
        size_t used = 0;
    };
    typedef std::unordered_map<InvariantVector, Info,
                               ImmutableVector_hash<size_t>,
                               ImmutableVector_equal<size_t>>
        Inv2Info;

    Invariants(Output &output, const BinaryFunction &table)
        : d_output(output), d_table(table) {}

    void calculate();

    Info &get(const InvariantVector &invv) {
        const auto i = d_invariants.find(invv);
        assert(d_invariants.end() != i);
        return i->second;
    }
    const Inv2Info &invariants() const { return d_invariants; }

  private:
    Output &d_output;
    const BinaryFunction &d_table;
    Inv2Info d_invariants;
};

/* Calculate how many steps are needed to do in the function in order to
 * loop around following the function fun, i.e. the oriented graph v ->
 * f(v).*/
class Looping {
  public:
    Looping(Output &output, const std::vector<size_t> &fun)
        : d_output(output), d_order(fun.size()), d_fun(fun),
          d_value(d_order, std::numeric_limits<std::size_t>::max()){};

    /* calculate looping for the element query_ix */
    size_t calc_loop(size_t query_ix);
    size_t operator()(size_t query_ix) { return calc_loop(query_ix); }

  private:
    Output &d_output;
    const size_t d_order;
    const std::vector<size_t> &d_fun;
    std::vector<size_t> d_value;
    bool has_val(size_t i) { return d_value[i] <= d_order; }
};

/* Calculate how many steps are needed to reach  a fixed element in the function
 * graph fun, i.e. the oriented graph v -> f(v).*/
class Distances {
  public:
    Distances(Output &output, const std::vector<size_t> &fun, size_t target)
        : d_output(output), d_order(fun.size()), d_fun(fun),
          d_infinity(d_order), d_undef(std::numeric_limits<std::size_t>::max()),
          d_distance(d_order, d_undef), d_target(target) {
        assert(d_infinity < d_undef);
        d_distance[d_target] = 0;
    }

    /* calculate distance for the element query_ix */
    size_t calc_distance(size_t query_ix);
    size_t operator()(size_t query_ix) { return calc_distance(query_ix); }

  private:
    Output &d_output;
    const size_t d_order;
    const std::vector<size_t> &d_fun;
    const size_t d_infinity;
    const size_t d_undef;
    std::vector<size_t> d_distance;
    size_t d_target;
    bool has_val(size_t i) { return d_distance[i] < d_undef; }
};

class BigInvariantCalculator {
  public:
    using InvariantType = BigInvariant::InvariantType;
    BigInvariantCalculator(Output &output, const BinaryFunction &table)
        : d_output(output), d_table(table) {}
    struct Info {
        std::set<size_t> elems;
    };
    typedef std::unordered_map<InvariantVector, Info,
                               ImmutableVector_hash<size_t>,
                               ImmutableVector_equal<size_t>>
        Inv2Info;

    void calculate();
    const Inv2Info &invariants() const { return d_invariants; }

  private:
    Output &d_output;
    const BinaryFunction &d_table;
    Inv2Info d_invariants;
    std::vector<size_t> d_occs;
    InvariantVector calculate(DiagInvariants &dgc, size_t i);
};

