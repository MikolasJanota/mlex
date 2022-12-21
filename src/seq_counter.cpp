/*
 * File:  SeqCounter.cc
 * Author:  mikolas
 * Created on:  Thu Dec 1 16:19:34 GMTST 2011
 * Copyright (C) 2011, Mikolas Janota
 */
#include "seq_counter.h"
#include "auxiliary.h"
/* #include "minisat_auxiliary.h" */
void SeqCounter::encode() {
    const size_t n = _literals.size();
    if (_tval == 0) {
        encode_all0();
        return;
    }
    if (_tval >= n)
        return;

    assert(n > 1);

    addClause(~v(1), s(1, 1));
    for (size_t j = 2; j <= _tval; ++j)
        addClause(~s(1, j));

    for (size_t i = 2; i < n; ++i) {
        addClause(~v(i), s(i, 1));
        addClause(~s(i - 1, 1), s(i, 1));

        for (size_t j = 2; j <= _tval; ++j) {
            addClause(~v(i), ~s(i - 1, j - 1), s(i, j));
            addClause(~s(i - 1, j), s(i, j));
        }
        addClause(~v(i), ~s(i - 1, _tval));
    }
    addClause(~v(n), ~s(n - 1, _tval));
}

void SeqCounter::encode_all0() {
    const size_t n = _literals.size();
    for (size_t i = 1; i <= n; ++i)
        addClause(~v(i));
}
