/*
 * File:  encoding_explicit.h
 * Author:  mikolas
 * Created on:  Tue Jul 18 16:43:33 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include "auxiliary.h" // for SATSPC
#include "binary_function.h"
#include "minisat/core/SolverTypes.h" // for Lit, operator~
#include "minisat/mtl/Vec.h"          // for vec
#include "minisat_ext.h"
#include "options.h"
#include "seq_counter.h"
#include <algorithm> // for max
#include <cstddef>   // for size_t
#include <memory>    // for allocator_traits<>::value_type
#include <sstream>
#include <string>        // for string, basic_string, hash
#include <tuple>         // for tuple
#include <unordered_map> // for unordered_map, operator!=, _No...
#include <utility>       // for pair
#include <vector>        // for vector

class EncodingExplicit {
  public:
    typedef std::tuple<size_t, size_t, size_t> Assignment;
    EncodingExplicit(Output &output, SATSPC::MiniSatExt &sat,
                     const BinaryFunction &table)
        : d_options(output.d_options), d_sat(sat), d_table(table) {
#ifndef NDEBUG
        LOGIPASIR(d_sat.set_representatives(&d_representatives););
#endif
        _encoding_pos_ls.capacity(4 * table.order());
        setup_perms();
        setup_vals();
    }

    void make_solution(BinaryFunction &table);
    void encode_bij();
    void encode_iso();
    SATSPC::Lit encode_less(const BinaryFunction &table);

  private:
    const Options &d_options;
    SATSPC::MiniSatExt &d_sat;
    const BinaryFunction &d_table;

    SATSPC::vec<SATSPC::Lit> _encoding_pos_ls;
    SATSPC::Lit encode_less_rec(const BinaryFunction &table, size_t r,
                                size_t c);

    inline void atm1(SATSPC::MiniSatExt &sat,
                     const std::vector<SATSPC::Lit> &literals) {
        LOGIPASIR(std::cout << "c atm1 { ";
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
        LOGIPASIR(std::cout << "c end atm1" << std::endl;);
    }

    inline void eq1(SATSPC::MiniSatExt &sat,
                    const std::vector<SATSPC::Lit> &literals) {
        atm1(sat, literals);
        sat.addClause(literals);
    }

    inline SATSPC::Lit lesslit() {
#ifdef NDEBUG
        return SATSPC::mkLit(d_sat.fresh());
#else
        static size_t counter = 0;
        std::stringstream sts;
        sts << "lt_" << counter++;
        return get_representative(sts.str());
#endif
    }

  public:
    inline SATSPC::Lit val(size_t r, size_t c, size_t v) {
        SATSPC::Lit &lit = d_vals[r][c][v];
        if (lit != SATSPC::lit_Undef)
            return lit;
#ifdef NDEBUG
        lit = SATSPC::mkLit(d_sat.fresh());
#else
        std::stringstream sts;
        sts << "[" << r << "," << c << "]=" << v;
        lit = get_representative(sts.str());
#endif
        return lit;
    }

    inline SATSPC::Lit perm(size_t dom, size_t rng) {
        auto &row = d_perms[dom];
        auto lit = row[rng];
        if (lit != SATSPC::lit_Undef)
            return lit;
#ifdef NDEBUG
        row[rng] = SATSPC::mkLit(d_sat.fresh());
#else
        std::stringstream sts;
        sts << dom << "->" << rng;
        row[rng] = get_representative(sts.str());
#endif
        return row[rng];
    }

    inline void setup_vals() {
        const auto n = d_table.order();
        d_vals.resize(n);
        for (size_t r = 0; r < n; r++) {
            d_vals[r].resize(n);
            for (size_t c = 0; c < n; c++)
                d_vals[r][c].resize(n, SATSPC::lit_Undef);
        }
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
    std::vector<std::vector<std::vector<SATSPC::Lit>>> d_vals;
};
