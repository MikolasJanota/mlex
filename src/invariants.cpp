/*
 * File:  invariants.cpp
 * Author:  mikolas
 * Created on:  Mon Jan 9 15:50:15 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "invariants.h"

void Invariants::calculate() {
    const auto n = d_table.order();
    InvariantCalculator calc(n);

    for (auto row = n; row--;) {
        calc.set_row(row);
        for (auto col = n; col--;)
            calc.set_val(col, d_table.get(row, col));
        InvariantVector iinv = calc.make_ivec();
        auto [it, _] = d_invariants.insert({iinv, Info()});
        it->second.original_rows.push_back(row);
    }
    if (d_options.verbose > 2) {
        for (const auto &[inv, info] : d_invariants) {
            d_output.comment(3) << "inv " << inv << " {";
            for (const auto k : info.original_rows)
                d_output.ccomment(3) << " " << k;
            d_output.ccomment(3) << " }" << std::endl;
        }
    }
}
