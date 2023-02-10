/*
 * File:  graph.h
 * Author:  mikolas
 * Created on:  Fri Feb 10 10:18:24 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "binary_function.h"
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

class Graph {
  public:
    Graph(const BinaryFunction &table) : d_table(table), d_id(0) {}
    void make();
    std::ostream &print_dot(std::ostream &out);
    std::ostream &print_nauty(std::ostream &out);

  private:
    const BinaryFunction &d_table;
    size_t d_id;
    std::vector<size_t> d_rowns, d_colns, d_valns;
    std::vector<std::vector<size_t>> d_cellns;
    std::vector<std::vector<size_t>> d_edges;

    void e(size_t a, size_t b) {
        assert(a < d_edges.size());
        d_edges[a].push_back(b);
    }
    void make_nodes();
    void make_edges();
    void make_colors();
    std::ostream &print_dot_node(std::ostream &out, size_t n, std::string label,
                                 const char *col);
};
