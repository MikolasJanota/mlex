/*
 * File:  lexmin_solver_explicit.h
 * Author:  mikolas
 * Created on:  Tue Jul 18 16:40:04 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "auxiliary.h" // for SATSPC
#include "encoding_explicit.h"
#include "lexmin_solver_base.h"
class LexminSolverExplicit : public LexminSolverBase {
  public:
    LexminSolverExplicit(Output &output, const BinaryFunction &table)
        : LexminSolverBase(output, table) {}

    virtual ~LexminSolverExplicit() {
        if (d_solution)
            delete d_solution;
    }

    virtual void solve() override;

    /* Make solution as compact BinaryFunction */
    virtual BinaryFunction *make_solution() override {
        auto *solution = d_solution;
        d_solution = nullptr;
        return solution;
    }

    /* Make solution as compact fun representation */
    virtual CompFunction make_solution_comp() override;

    virtual void set_diag(const std::vector<size_t> &) override {
        assert(false);
    }

  private:
    std::unique_ptr<SATSPC::MiniSatExt> d_sat;
    std::unique_ptr<EncodingExplicit> d_encoding;
    BinaryFunction *d_solution = nullptr;
};
