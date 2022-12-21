/*
 * File:  BinaryFunction.cpp
 * Author:  mikolas
 * Created on:  Mon Aug 10 13:54:06 WEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#include<iostream>
#include"binary_function.h"

void BinaryFunction::print(std::ostream& out) const {
    for (size_t i = 0; i < order(); i++) {
        for (size_t j = 0; j < order(); j++) {
            out << get(i, j);
            if (i + 1 < order() || j + 1 < order())
                out << ",";
        }
        out << std::endl;
    }
}
