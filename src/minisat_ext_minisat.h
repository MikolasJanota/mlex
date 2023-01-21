/*
 * File:   MiniSatExt.hh
 * Author: mikolas
 *
 * Created on November 29, 2010, 5:40 PM
 */
#pragma once
#include "auxiliary.h"
/* #define USE_MINISATSIMP */
#ifdef USE_MINISATSIMP
#include "minisat/simp/SimpSolver.h"
#else
#include "minisat/core/Solver.h"
#endif
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
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

#ifdef USE_MINISATSIMP
    inline void setFrozen(Var v, bool b) { _solver.setFrozen(v, b); }
#else
    inline void setFrozen(Var, bool) {}
#endif

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
        int i = 0;
        for (const auto l : cl)
            ls[i++] = l;
        return _solver.addClause_(ls);
    }

    inline Minisat::lbool get_model_value(Minisat::Var v) const {
        return v < _solver.model.size() ? _solver.model[v] : l_Undef;
    }

    inline Minisat::lbool eval_lit(const Minisat::Lit &l) const {
        const Minisat::lbool lval = get_model_value(var(l));
        return lval == Minisat::l_Undef
                   ? Minisat::l_Undef
                   : (Minisat::sign(l) == (lval == Minisat::l_False)
                          ? Minisat::l_True
                          : Minisat::l_False);
    }

  private:
#ifdef USE_MINISATSIMP
    SimpSolver _solver;
#else
    Solver _solver;
#endif
};

inline void MiniSatExt::new_variables(Var max_id) {
    const int target_number = (int)max_id + 1;
    while (_solver.nVars() < target_number)
        _solver.newVar();
}

} // namespace SATSPC
