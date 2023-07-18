/*
 * File:  lexmin_solver_base.h
 * Author:  mikolas
 * Created on:  Tue Jul 18 14:37:43 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "binary_function.h" // for BinaryFunction
#include "comp_function.h"
#include "options.h"
#include <vector>

class StatisticsManager;

class LexminSolverBase {
  public:
    LexminSolverBase(Output &output, const BinaryFunction &table)
        : d_output(output), d_options(output.d_options),
          d_statistics(output.d_statistics), d_table(table) {}

    virtual ~LexminSolverBase() {}
    virtual void solve() = 0;

    /* Make solution as compact BinaryFunction */
    virtual BinaryFunction *make_solution() = 0;

    /* Make solution as compact fun representation */
    virtual CompFunction make_solution_comp() = 0;
    virtual void set_diag(const std::vector<size_t> &diag) = 0;

  protected:
    Output &d_output;
    const Options &d_options;
    StatisticsManager &d_statistics;
    const BinaryFunction &d_table;

    bool d_is_solved = false;

    inline std::ostream &comment(int level = 0) {
        return d_output.comment(level);
    }

    inline std::ostream &ccomment(int level = 0) {
        return d_output.ccomment(level);
    }
};
