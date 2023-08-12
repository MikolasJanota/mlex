/*
 * File:  color_invariants.cpp
 * Author:  mikolas
 * Created on:  Fri Aug 11 13:36:40 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */

#include "color_invariants.h"
#include "union_find.h"
#include <cstddef>
#include <map>
#include <set>

void ColorInvariantCalculator::print_node_invariant(
    size_t node, const std::vector<size_t> &inv) const {
    assert(inv.size() == d_node_invariant_size);
    assert(inv.size() == d_color_count + 3);
    d_output.comment(3) << d_row << "." << node << ":[";
    for (size_t i = 0; i < d_color_count; ++i)
        d_output.ccomment(3) << (i ? " " : "") << "C" << i << ":" << inv[i];
    d_output.ccomment(3) << "|NC" << inv[d_color_count];
    d_output.ccomment(3) << "|MC" << inv[d_color_count + 1];
    d_output.ccomment(3) << "|" << inv[d_color_count + 2];
    d_output.ccomment(3) << "]" << std::endl;
}

void ColorInvariantCalculator::make_invariants(
    std::vector<std::vector<size_t>> &invariants) {
    invariants.clear();
    invariants.resize(d_n, std::vector<size_t>(d_node_invariant_size, 0));
    for (size_t i = 0; i < d_n; ++i) {
        const auto next = d_values[i];
        assert(next < d_n);
        const auto curr_color = d_colors[i];
        const auto next_color = d_colors[next];
        assert(curr_color < d_color_count);
        assert(next_color < d_color_count);
        invariants[next][curr_color]++;
        invariants[i][d_color_count] = next_color;
        invariants[i][d_color_count + 1] = curr_color;
        invariants[i][d_color_count + 2] = (i == d_row) ? 1 : 0;
    }
}

void ColorInvariantCalculator::make_map(InvMap &inv_map) {
    std::vector<std::vector<size_t>> invariants;
    make_invariants(invariants);
    for (size_t i = 0; i < d_n; ++i) {
        const auto inv = InvariantVector(invariants[i]);
        inv_map[inv].insert(i);
    }
}

void ColorInvariantCalculator::calculate() {
    std::vector<std::vector<size_t>> invariants;
    make_invariants(invariants);
    print_vec(d_output.comment(3) << d_row << ":", d_values) << std::endl;
    d_invariants.resize(d_n);
    for (size_t i = 0; i < d_n; ++i) {
        const auto inv = InvariantVector(invariants[i]);
        d_invariants[i] = inv;
        d_frequencies[inv] = 0;
    }
    if (d_output.d_options.verbose > 2)
        for (size_t i = 0; i < d_n; ++i)
            print_node_invariant(i, invariants[i]);
    for (const auto &inv : d_invariants)
        (d_frequencies[inv])++;
}

void ColorInvariantManager::build_colors() {
    UnionFind<size_t> uf;
    const auto n = d_table.order();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            if (d_allowed_mapping.allowed(i, j)) {
                d_output.comment(3) << "uf:" << i << "->" << j << std::endl;
                uf.unify(to_int(i, true), to_int(j, false));
            }

    std::set<size_t> representatives;
    for (size_t i = 0; i < n; ++i)
        representatives.insert(uf.get_repr(to_int(i, true)));
    std::map<size_t, size_t> repr2color;
    d_color_count = 0;
    d_colors_src.resize(n, -1);
    d_colors_dst.resize(n, -1);
    for (const auto r : representatives)
        repr2color[r] = d_color_count++;

    for (size_t i = 0; i < n; ++i) {
        const auto r = uf.get_repr(to_int(i, true));
        d_colors_src[i] = repr2color.at(r);
        d_output.comment(3)
            << "colsrc:" << i << "->" << d_colors_src[i] << std::endl;
    }
    for (size_t i = 0; i < n; ++i) {
        const auto r = uf.get_repr(to_int(i, false));
        d_colors_dst[i] = repr2color.at(r);
        d_output.comment(3)
            << "coldst:" << i << "->" << d_colors_dst[i] << std::endl;
    }
}

bool ColorInvariantManager::add_row(size_t dst_row,
                                    const BinaryFunction &table) {
    const auto n = d_table.order();
    assert(n == table.order());
    d_colors_src.clear();
    d_colors_dst.clear();
    d_src_row_color_invariants.clear();
    d_dst_row_color_invariants.clear();
    d_row_inv_src.clear();
    d_row_inv_dst.clear();
    if (d_colors_dst.empty())
        build_colors();

    ColorInvariantCalculator srcc(d_output, d_color_count, d_colors_src);
    d_src_row_color_invariants.clear();
    d_row_inv_src.clear();
    for (size_t row = 0; row < n; ++row) {
        srcc.set_row(row);
        for (size_t col = 0; col < n; ++col)
            srcc.set_val(col, d_table.get(row, col));
        srcc.calculate();
        const auto inv = srcc.inv();
        d_output.comment(3) << "src row inv:" << row << std::endl;
        print(inv);
        d_src_row_color_invariants[inv].insert(row);
        d_row_inv_src.push_back(inv);
    }

    ColorInvariantCalculator dstc(d_output, d_color_count, d_colors_dst);
    d_dst_row_color_invariants.clear();
    d_row_inv_dst.clear();
    for (size_t row = 0; row <= dst_row; ++row) {
        dstc.set_row(row);
        for (size_t col = 0; col < n; ++col)
            dstc.set_val(col, table.get(row, col));
        dstc.calculate();
        const auto inv = dstc.inv();
        d_output.comment(3) << "dst row inv:" << row << std::endl;
        print(inv);
        d_dst_row_color_invariants[inv].insert(row);
        d_row_inv_dst.push_back(inv);
    }
    return true;
}

void ColorInvariantManager::print(const Invariant &inv) {
    for (const auto &[ninv, count] : inv)
        d_output.comment(3) << ninv << ":" << count << std::endl;
}

