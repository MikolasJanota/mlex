/*
 * File:  comp_function.cpp
 * Author:  mikolas
 * Created on:  Fri Jan 27 15:41:20 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#include "comp_function.h"
#include <cstdlib>
std::ostream &CompFunction::print_gap(std::ostream &output) const {
    if (d_arity != 2) {
        std::cerr << "arity other than 2 not handled yet" << std::endl;
        exit(1);
    }

    const auto n = order();
    output << "[" << std::endl;
    CompFunctionReader r(*this);
    const auto sz = calculate_cell_count(n, d_arity);
    for (size_t i = 0; i < sz; i++) {
        const auto row = i / n;
        const auto col = i % n;
        if (col == 0)
            output << "   [";
        output << (1 + r.next());
        if ((col + 1) == n)
            output << "]" << (row + 1 == n ? "" : ",") << std::endl;
        else
            output << ",";
    }
    return output << "]";
}

std::ostream &CompFunction::print_mace(std::ostream &output,
                                       const std::string &info) const {
    if (d_arity != 2) {
        std::cerr << "arity other than 2 not handled yet" << std::endl;
        exit(1);
    }

    const auto n = order();
    output << "interpretation( " << n << ", [" << info << "], [" << std::endl;
    output << "  function(" << get_name() << "(_,_), [" << std::endl;
    CompFunctionReader r(*this);
    const auto sz = calculate_cell_count(n, d_arity);
    for (size_t i = 0; i < sz; i++) {
        /* const auto row = i / n; */
        const auto col = i % n;
        if (col == 0)
            output << "    ";
        output << r.next();
        if (i + 1 < sz) {
            output << ",";
            if ((col + 1) == n)
                output << std::endl;
        }
    }
    return output << "])]).";
}

