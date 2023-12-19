/*
 * File:  BinaryFunction.cpp
 * Author:  mikolas
 * Created on:  Mon Aug 10 13:54:06 WEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#include "binary_function.h"
#include <iomanip> // std::setw
#include <iostream>

std::ostream &BinaryFunction::print(std::ostream &out, const char *pre) const {
    out << pre;
    for (size_t j = 0; j < order(); j++)
        out << ((j == 0) ? "  " : " ") << j;
    out << std::endl;

    for (size_t i = 0; i < order(); i++) {
        out << pre << i << ":";
        for (size_t j = 0; j < order(); j++) {
            out << get(i, j);
            if (i + 1 < order() || j + 1 < order())
                out << ",";
        }
        out << std::endl;
    }
    return out;
}

std::ostream &BinaryFunction::print_gap(std::ostream &output) const {
    const auto n = order();
    output << "[ " << std::endl;
    for (size_t row = 0; row < n; row++) {
        output << "[";
        for (size_t col = 0; col < n; col++)
            output << (col ? "," : "") << (get(row, col) + 1);
        output << "]" << ((row + 1) < n ? "," : "") << std::endl;
    }
    return output << "]";
}

std::ostream &BinaryFunction::print_mace(std::ostream &output) const {
    const auto n = order();
    output << "interpretation( " << n << ", [" << get_additional_info()
           << "], [" << std::endl;
    output << "  function(" << get_name() << "(_,_), [" << std::endl;
    for (size_t row = 0; row < n; row++) {
        output << "    ";
        for (size_t col = 0; col < n; col++) {
            if (col)
                output << ",";
            if (is_set(row, col))
                output << get(row, col);
            else
                output << -1;
        }
        if (row + 1 < n)
            output << "," << std::endl;
    }
    return output << "])]).";
}
