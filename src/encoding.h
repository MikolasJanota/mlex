/*
 * File:  encoding.h
 * Author:  mikolas
 * Created on:  Tue Dec 13 12:20:30 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include "binary_function.h"
#include "minisat/core/SolverTypes.h" // for Lit, operator~
#include "minisat/mtl/Vec.h"          // for vec
#include "minisat_ext.h"
#include "options.h"
#include "seq_counter.h"
#include "auxiliary.h" // for SATSPC
#include <memory>      // for allocator_traits<>::value_type
#include <sstream>
#include <cstddef>      // for size_t
#include <string>        // for string, basic_string, hash
#include <tuple>         // for tuple
#include <unordered_map> // for unordered_map, operator!=, _No...
#include <utility>       // for pair
#include <vector>        // for vector

class Encoding {
  public:
    typedef std::tuple<size_t, size_t, size_t> Assignment;
    Encoding(Output &output, SATSPC::MiniSatExt &sat,
             const BinaryFunction &table)
        : d_options(output.d_options), d_sat(sat), d_table(table) {
#ifndef NDEBUG
        LOGIPASIR(d_sat.set_representatives(&d_representatives););
#endif
        setup_perms();
    }

    void print_solution(std::ostream &output);

    void encode_bij();
    void encode(const std::vector<Assignment> &assignments);
    void encode_pos(const Assignment &val, SATSPC::Lit selector);
    void encode_shot(const std::pair<size_t, size_t> &cell,
                     const std::vector<size_t> &vals, SATSPC::Lit selector);

  private:
    const Options &d_options;
    SATSPC::MiniSatExt &d_sat;
    const BinaryFunction &d_table;

    SATSPC::vec<SATSPC::Lit> _encoding_pos_ls;

    inline void atm1(SATSPC::MiniSatExt &sat,
                     const std::vector<SATSPC::Lit> &literals) {
        LOGIPASIR(std::cout << "atm1 { ";
                  for (const auto &literal
                       : literals) sat.print_literal(std::cout, literal)
                  << " ";
                  std::cout << "}" << std::endl;);
        if (literals.size() < d_options.seq_counter_lits) {
            for (size_t i = 0; (i + 1) < literals.size(); i++)
                for (size_t j = i + 1; j < literals.size(); j++)
                    sat.addClause(~literals[i], ~literals[j]);
        } else {
            SeqCounter(sat, literals, 1).encode();
        }
        LOGIPASIR(std::cout << "end atm1" << std::endl;);
    }

    inline void eq1(SATSPC::MiniSatExt &sat,
                    const std::vector<SATSPC::Lit> &literals) {
        atm1(sat, literals);
        sat.addClause(literals);
    }

  public:
    inline SATSPC::Lit perm(size_t dom, size_t rng) {
        auto &row = d_perms[dom];
        auto lit = row[rng];
        if (lit != SATSPC::lit_Undef)
            return lit;
#ifdef NDEBUG
        row[rng] = SATSPC::mkLit(d_sat.fresh());
#else
        std::stringstream sts;
        sts << "p_" << dom << "_" << rng;
        row[rng] = get_representative(sts.str());
#endif
        return row[rng];
    }

    inline void setup_perms() {
        const auto n = d_table.order();
        d_perms.resize(n);
        for (size_t dom = 0; dom < n; dom++)
            d_perms[dom].resize(n, SATSPC::lit_Undef);
    }

#ifndef NDEBUG
    std::unordered_map<std::string, SATSPC::Lit> d_representatives;

    SATSPC::Lit get_representative(const std::string &name) {
        const auto index = d_representatives.find(name);
        if (index != d_representatives.end())
            return index->second;
        const auto representative = SATSPC::mkLit(d_sat.fresh());
        d_representatives.insert(index, {name, representative});
        LOGIPASIR(d_sat.d_inverse_representatives->insert(
            {SATSPC::var(representative), name}););
        return representative;
    }
#endif

  private:
    std::vector<std::vector<SATSPC::Lit>> d_perms;
};
