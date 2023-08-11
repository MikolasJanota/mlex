/*
 * File:  color_invariants.h
 * Author:  mikolas
 * Created on:  Fri Aug 11 13:19:00 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "allowed_mapping.h"
#include "binary_function.h"
#include "immutable_vector.h"
#include "invariants.h"
#include <cstddef>
#include <list>
#include <ostream>
#include <unordered_map>
#include <vector>

/* For a given colored row calculate invariants. */
class ColorInvariantCalculator {
  public:
    typedef std::unordered_map<InvariantVector, size_t,
                               ImmutableVector_hash<size_t>,
                               ImmutableVector_equal<size_t>>
        InvHistogram;

    ColorInvariantCalculator(size_t color_count,
                             const std::vector<size_t> &colors)
        : d_n(colors.size()), d_color_count(color_count),
          d_node_invariant_size(d_color_count + 2), d_colors(colors),
          d_values(d_n, -1) {
        assert(d_colors.size() == d_n);
    }

    /* resets the class for constructing invariant for a new row */
    void set_row(size_t row) {
        d_row = row;
        d_values.clear();
        d_values.resize(d_n, -1);
        d_invariants.clear();
        d_frequencies.clear();
    }

    /* register a value val in column col */
    void set_val(size_t col, size_t val) { d_values[col] = val; }

    /* process current row */
    void calculate();

    InvHistogram inv() { return d_frequencies; };

  protected:
    const size_t d_n;

    // number of colors we are currently using
    const size_t d_color_count;

    // size of invariants vector for each node
    const size_t d_node_invariant_size;

    // coloring of elements
    const std::vector<size_t> &d_colors;

    // values of elements in the row
    std::vector<size_t> d_values;

    // invariant for each element
    std::vector<InvariantVector> d_invariants;

    // node invariant frequencies
    InvHistogram d_frequencies;

    // number of current row
    size_t d_row;
};

class InvHistogram_eq {
  public:
    bool operator()(const ColorInvariantCalculator::InvHistogram &h1,
                    const ColorInvariantCalculator::InvHistogram &h2) const {
        if (h1.size() != h2.size())
            return false;
        for (const auto &[inv, count] : h1) {
            const auto i = h2.find(inv);
            if (i == h2.end() || i->second != count)
                return false;
        }
        return true;
    }
};
class InvHistogram_hash {
  public:
    size_t operator()(const ColorInvariantCalculator::InvHistogram &h) const {
        size_t rv = h.size();
        for (const auto &[inv, count] : h)
            rv ^= inv.hash_code() ^ count;
        return rv;
    }
};

struct ColorInvariantManager {
    explicit ColorInvariantManager(Output &output, const BinaryFunction &f)
        : d_output(output), d_table(f), d_allowed_mapping(d_table.order()) {}
    typedef ColorInvariantCalculator::InvHistogram Invariant;
    typedef std::unordered_map<Invariant, std::set<size_t>, InvHistogram_hash,
                               InvHistogram_eq>
        InvMap;

    Output &d_output;
    const BinaryFunction &d_table;

    size_t d_color_count = 0;
    ;
    std::vector<size_t> d_colors_src;
    std::vector<size_t> d_colors_dst;
    InvMap d_src_row_color_invariants;
    InvMap d_dst_row_color_invariants;
    AllowedMapping d_allowed_mapping;

    bool add_row(size_t row, const BinaryFunction &table);
    void build_colors();
    static size_t to_int(size_t node, bool out) {
        return (2 * node) + (out ? 1 : 0);
    }

    void print(const Invariant &inv);
};
