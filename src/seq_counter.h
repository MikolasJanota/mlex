/*
 * File:  SeqCounter.hh
 * Author:  mikolas
 * Created on:  Thu Dec 1 15:53:16 GMTST 2011
 * Copyright (C) 2011, Mikolas Janota
 */
#pragma once
#define SEQCOUNTER_DBG(t)
#include "auxiliary.h"
#include "minisat_ext.h"
#include <iostream>
#include <vector>

class SeqCounter {
  public:
    SeqCounter(SATSPC::MiniSatExt &solver,
               const std::vector<SATSPC::Lit> &pliterals, size_t ptval)
        : _solver(solver), _literals(pliterals), _tval(ptval) {
        if (!_literals.empty()) {
            auto needed = (_literals.size() - 1) * _tval;
            while (needed--)
                aux.push_back(SATSPC::mkLit(_solver.fresh()));
        }
    }

    void encode();
    void encode_all0();

  private:
    SATSPC::MiniSatExt &_solver;
    const std::vector<SATSPC::Lit> &_literals;
    const size_t _tval;
    std::vector<SATSPC::Lit> aux;

    inline SATSPC::Lit s(size_t v_index, size_t j) {
        assert(v_index > 0);
        assert(v_index < _literals.size());
        assert(j >= 1);
        assert(j <= _tval);
        const size_t aux_index = (v_index - 1) * _tval + j - 1;
        assert(aux_index < aux.size());
        return aux[aux_index];
    }

    inline SATSPC::Lit v(size_t v_index) {
        assert(v_index > 0);
        assert(v_index <= _literals.size());
        return _literals[v_index - 1];
    }

    inline void addClause(SATSPC::Lit p);
    inline void addClause(SATSPC::Lit p, SATSPC::Lit q);
    inline void addClause(SATSPC::Lit p, SATSPC::Lit q, SATSPC::Lit r);
};

inline void SeqCounter::addClause(SATSPC::Lit p) { _solver.addClause(p); }

inline void SeqCounter::addClause(SATSPC::Lit p, SATSPC::Lit q) {
    _solver.addClause(p, q);
}

inline void SeqCounter::addClause(SATSPC::Lit p, SATSPC::Lit q, SATSPC::Lit r) {
    _solver.addClause(p, q, r);
}
