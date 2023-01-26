/*
 * File:  BinaryFunction.cpp
 * Author:  mikolas
 * Created on:  Mon Aug 10 13:54:06 WEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#include "binary_function.h"
#include <climits>
#include <iostream>

void BinaryFunction::print(std::ostream &out) const {
    for (size_t i = 0; i < order(); i++) {
        for (size_t j = 0; j < order(); j++) {
            out << get(i, j);
            if (i + 1 < order() || j + 1 < order())
                out << ",";
        }
        out << std::endl;
    }
}
void BinaryFunction::print_gap(std::ostream &output) {
    const auto n = order();
    output << "[ " << std::endl;
    for (size_t row = 0; row < n; row++) {
        output << "[ ";
        for (size_t col = 0; col < n; col++) {
            if (col)
                output << " , ";
            output << get(row, col) + 1;
        }
        output << " ]," << std::endl;
        ;
    }
    output << "]" << std::endl;
}

void CompFunction::print_mace(std::ostream &output) const {
    if (d_arity != 2) {
        std::cerr << "arity other than 2 not handled yet" << std::endl;
        exit(1);
    }

    const auto n = order();
    output << "interpretation( " << n << ", [" << get_additional_info()
           << "], [" << std::endl;
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
    output << "])])." << std::endl;
}

void BinaryFunction::print_mace(std::ostream &output) {
    const auto n = order();
    output << "interpretation( " << n << ", [" << get_additional_info()
           << "], [" << std::endl;
    output << "  function(" << get_name() << "(_,_), [" << std::endl;
    for (size_t row = 0; row < n; row++) {
        output << "    ";
        for (size_t col = 0; col < n; col++) {
            if (col)
                output << ",";
            output << get(row, col);
        }
        if (row + 1 < n)
            output << "," << std::endl;
        ;
    }
    output << "])])." << std::endl;
}

