/*
 * File:  graph.cpp
 * Author:  mikolas
 * Created on:  Fri Feb 10 10:18:33 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "graph.h"
#include "auxiliary.h"
#include <cstddef>
#include <string>
#include <vector>

void Graph::make_nodes() {
    const auto n = d_table.order();
    for (size_t _ = 0; _ < n; ++_)
        d_rowns.push_back(d_id++);
    for (size_t _ = 0; _ < n; ++_)
        d_colns.push_back(d_id++);
    d_cellns.resize(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t _ = 0; _ < n; ++_) {
            d_cellns[i].push_back(d_id++);
        }
    }
    for (size_t _ = 0; _ < n; ++_)
        d_valns.push_back(d_id++);
}

std::ostream &print_vec(std::ostream &out, const std::vector<size_t> &vs) {
    bool f = true;
    for (const auto i : vs) {
        if (!f)
            out << ",";
        out << i;
        f = false;
    }
    return out;
}

std::ostream &Graph::print_nauty(std::ostream &out) {
    out << "l=0 -m -a n=" << d_id << " g\n";
    const auto esz = d_edges.size();
    bool fstln = true;
    for (size_t i = 0; i < esz; i++) {
        const auto &ns = d_edges[i];
        if (ns.empty())
            continue;
        if (!fstln)
            out << ",\n";
        fstln = false;
        out << i << ":";
        for (auto v : ns)
            out << " " << v;
    }
    out << ".\n";
    out << "f=[ ";
    print_vec(out, d_rowns) << " | ";
    print_vec(out, d_colns) << " | ";
    print_vec(out, d_valns) << " | ";

    bool fst = true;
    for (const auto &vs : d_cellns) {
        for (const auto i : vs) {
            if (!fst)
                out << ",";
            out << i;
            fst = false;
        }
    }
    out << " ]\n";

    return out << "c x z b q" << std::endl;
}

std::ostream &Graph::print_dot_node(std::ostream &out, size_t n,
                                    std::string label, const char *col) {
    return out << "    " << n << " [color=" << col << ",label=\"" << label
               << "\",style=filled];\n";
}

std::ostream &Graph::print_dot(std::ostream &out) {
    const auto n = d_table.order();
    out << "graph G {\n";
    for (size_t row = 0; row < n; ++row)
        print_dot_node(out, d_rowns[row], "r" + std::to_string(row), "red");
    for (size_t col = 0; col < n; ++col)
        print_dot_node(out, d_colns[col], "c" + std::to_string(col), "green");
    for (size_t val = 0; val < n; ++val)
        print_dot_node(out, d_valns[val], std::to_string(val), "steelblue3");
    for (size_t r = 0; r < n; ++r)
        for (size_t c = 0; c < n; ++c)
            print_dot_node(out, d_cellns[r][c],
                           "a" + std::to_string(r) + "_" + std::to_string(c),
                           "goldenrod2");
    for (size_t i = 0; i < d_id; ++i) {
        if (d_edges[i].empty())
            continue;
        out << "    " << i << " -- {";
        for (auto j : d_edges[i])
            out << " " << j;
        out << " };\n";
    }
    return out << "}";
}

void Graph::make_edges() {
    const auto n = d_table.order();
    d_edges.resize(d_id);
    for (size_t row = 0; row < n; ++row)
        for (size_t col = 0; col < n; ++col) {
            const auto valn = d_valns[d_table.get(row, col)];
            const auto rown = d_rowns[row];
            const auto coln = d_colns[col];
            const auto celln = d_cellns[row][col];
            e(celln, coln);
            e(celln, rown);
            e(celln, valn);
        }
    for (size_t i = 0; i < n; ++i) {
        const auto valn = d_valns[i];
        const auto rown = d_rowns[i];
        const auto coln = d_colns[i];
        e(valn, rown);
        e(valn, coln);
        e(coln, rown);
    }
}

void Graph::make() {
    make_nodes();
    make_edges();
}
