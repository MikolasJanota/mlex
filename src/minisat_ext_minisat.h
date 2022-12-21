/*
 * File:   MiniSatExt.hh
 * Author: mikolas
 *
 * Created on November 29, 2010, 5:40 PM
 */
#pragma once
#include "auxiliary.h"
/* #include "minisat/core/Solver.h" */
#include "minisat/simp/SimpSolver.h"
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#define SATSPC Minisat
namespace SATSPC {
std::ostream &operator<<(std::ostream &outs, Lit lit);

class MiniSatExt {
  public:
    inline const Minisat::LSet &conflict() { return _solver.conflict; }
    inline const Minisat::vec<Minisat::lbool> &model() { return _solver.model; }
    inline Var fresh() { return _solver.newVar(); }
    /* inline void bump(Var var) { _solver.varBumpActivity(var); } */
    inline void new_variables(Var max_id);

    inline bool solve() { return _solver.solve(); }
    inline bool solve(const vec<Lit> &assumps) {
        return _solver.solve(assumps);
    }

    void freezeVar(Var v) { _solver.freezeVar(v); }
    void thaw() { _solver.thaw(); }

    inline bool addClause(Lit a) { return _solver.addClause(a); }
    inline bool addClause(Lit a, Lit b) { return _solver.addClause(a, b); }
    inline bool addClause(Lit a, Lit b, Lit c) {
        return _solver.addClause(a, b, c);
    }

    inline bool addClause_(SATSPC::vec<SATSPC::Lit> &cl) {
        return _solver.addClause_(cl);
    }

    inline bool addClause(const std::vector<SATSPC::Lit> &cl) {
        SATSPC::vec<SATSPC::Lit> ls(cl.size());
        int                      i = 0;
        for (const auto l : cl)
            ls[i++] = l;
        return _solver.addClause_(ls);
    }

    inline Minisat::lbool get_model_value(Minisat::Var v) const {
        return v < _solver.model.size() ? _solver.model[v] : l_Undef;
    }

  private:
    SimpSolver _solver;
};

inline void MiniSatExt::new_variables(Var max_id) {
    const int target_number = (int)max_id + 1;
    while (_solver.nVars() < target_number)
        _solver.newVar();
}

} // namespace SATSPC
