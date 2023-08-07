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
        if (d_sol)
            delete d_sol;
    }

    virtual void solve() override;

    /* Make solution as BinaryFunction */
    virtual BinaryFunction *make_solution() override {
        auto *solution = d_sol;
        d_sol = nullptr;
        return solution;
    }

    /* Make solution as compact fun representation */
    virtual CompFunction make_solution_comp() override;

    virtual void set_diag(const std::vector<size_t> &) override {
        assert(false);
    }

    void make_permutation(std::vector<size_t> &perm);

  private:
    std::unique_ptr<SATSPC::MiniSatExt> d_sat;
    std::unique_ptr<EncodingExplicit> d_encoding;
    BinaryFunction *d_sol = nullptr;
    bool run_sat(Minisat::vec<Minisat::Lit> &assumps);
    bool run_sat();
    void opt1stRow();
    void solve_lowering();
    void solve_binary();
    void solve_left_to_right();
};
