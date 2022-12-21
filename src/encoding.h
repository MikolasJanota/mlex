/*
 * File:  encoding.h
 * Author:  mikolas
 * Created on:  Tue Dec 13 12:20:30 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include "binary_function.h"
#include "minisat_ext.h"
#include "options.h"
#include "seq_counter.h"
#include <sstream>

inline void atm1(SATSPC::MiniSatExt &            sat,
                 const std::vector<SATSPC::Lit> &literals) {
    LOGIPASIR(std::cout << "atm1 { ";
              for (const auto &literal
                   : literals) sat.print_literal(std::cout, literal)
              << " ";
              std::cout << "}" << std::endl;);
    if (literals.size() < 10) {
        for (size_t i = 0; (i + 1) < literals.size(); i++)
            for (size_t j = i + 1; j < literals.size(); j++)
                sat.addClause(~literals[i], ~literals[j]);
    } else {
        SeqCounter(sat, literals, 1).encode();
    }
    LOGIPASIR(std::cout << "end atm1" << std::endl;);
}

inline void eq1(SATSPC::MiniSatExt &            sat,
                const std::vector<SATSPC::Lit> &literals) {
    atm1(sat, literals);
    sat.addClause(literals);
}

class Encoding {
  public:
    typedef std::tuple<size_t, size_t, size_t> Assignment;
    Encoding(const Options &options, SATSPC::MiniSatExt &sat,
             const BinaryFunction &table)
        : d_options(options), d_sat(sat), d_table(table) {
        LOGIPASIR(d_sat.set_representatives(&d_representatives););
    }

    void print_solution(std::ostream &output);

    void encode_bij();
    void encode(const std::vector<Assignment> &assignments);
    void encode_pos(const Assignment &val, SATSPC::Lit selector);

  private:
    const Options &                              d_options;
    SATSPC::MiniSatExt &                         d_sat;
    const BinaryFunction &                       d_table;
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

  public:
    SATSPC::Lit perm(size_t a, size_t b) {
        std::stringstream sts;
        sts << "p_" << a << "_" << b;
        return get_representative(sts.str());
    }
};
