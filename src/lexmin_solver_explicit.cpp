/*
 * File:  lexmin_solver_explicit.cpp
 * Author:  mikolas
 * Created on:  Tue Jul 18 18:02:15 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "lexmin_solver_explicit.h"
#include "binary_function.h"
#include "encoding_explicit.h"
#include <cassert>
using SATSPC::Lit;

void LexminSolverExplicit::solve() {
    const auto n = d_table.order();
    d_sat = std::make_unique<SATSPC::MiniSatExt>();
    d_encoding = std::make_unique<EncodingExplicit>(d_output, *d_sat, d_table);
    d_encoding->encode_bij();
    d_encoding->encode_iso();
    if (d_solution)
        delete d_solution;
    d_solution = new BinaryFunction(n);
    d_solution->set(d_table);
    Lit last = d_encoding->encode_less(d_table);
    Minisat::vec<Minisat::Lit> assumps(1, last);
    while (d_sat->solve(assumps)) {
        d_sat->addClause(~last);
        d_encoding->make_solution(*d_solution);
        last = d_encoding->encode_less(d_table);
        assumps[0] = last;
    }
}

CompFunction LexminSolverExplicit::make_solution_comp() {
    const auto n = d_table.order();
    CompFunctionBuilder b(n, 2);
    for (size_t r = 0; r < n; r++)
        for (size_t c = 0; c < n; c++)
            b.push(d_solution->get(r, c));
    return b.make();
}
