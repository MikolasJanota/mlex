/*
 * File:  lit_funcs.h
 * Author:  mikolas
 * Created on:  Tue Nov 9 19:01:07 CET 2021
 * Copyright (C) 2021, Mikolas Janota
 */
#ifndef LIT_FUNCS_H_18168
#define LIT_FUNCS_H_18168
#include "minisat/core/SolverTypes.h"
class Lit_equal {
public:
  inline bool operator () (const Minisat::Lit& l1,const Minisat::Lit& l2) const { return l1==l2; }
};

class Lit_hash {
public:
  inline size_t operator () (const Minisat::Lit& l) const { return Minisat::toInt(l); }
};


#endif /* LIT_FUNCS_H_18168 */
