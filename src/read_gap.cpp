/*
 * File:  ReadGAP.cpp
 * Author:  mikolas
 * Created on:  Fri Oct 2 09:31:24 CEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#include "read_gap.h"
#include "auxiliary.h"
#include "fmtutils.hh"
#include <cstdlib>
#include <iostream>
#include <vector>
using std::vector;

ReadGAP::ReadGAP(gzFile &input_file) : _input_file(input_file) {}

static void match_char(StreamBuffer &sb, char c) {
    skipWhitespace(sb);
    while (*sb == '#') {
        skipLine(sb);
        skipWhitespace(sb);
    }
    if (*sb != c) {
        std::cerr << c << " expected" << std::endl;
        exit(EXIT_FAILURE);
    }
    ++sb;
}

void ReadGAP::read() {
    StreamBuffer in(_input_file);
    match_char(in, '[');

    vector<int>         row;
    vector<vector<int>> rows;
    size_t              order = -1;

    while (true) {
        skipWhitespace(in);
        match_char(in, '[');
        row.clear();
        while (true) {
            skipWhitespace(in);
            row.push_back(parseInt(in));
            skipWhitespace(in);
            if (*in == ']') {
                ++in;
                break;
            }
            match_char(in, ',');
        }
        if (rows.size() && row.size() != order) {
            std::cerr << "invalid row size " << std::endl;
            exit(EXIT_FAILURE);
        }
        order = row.size();
        rows.push_back(row);
        skipWhitespace(in);
        if (*in == ']') {
            ++in;
            break;
        }
        match_char(in, ',');
    }
    if (rows.size() != order) {
        std::cerr << "invalid number of rows" << std::endl;
        exit(EXIT_FAILURE);
    }
    _f = std::make_unique<BinaryFunction>(order);
    for (size_t i = 0; i < order; i++) {
        for (size_t j = 0; j < order; j++) {
            const auto number = rows[i][j];
            if (number <= 0 || number > static_cast<int>(order)) {
                std::cerr << "invalid number " << number << std::endl;
                exit(EXIT_FAILURE);
            }
            _f->set(i, j, number - 1);
        }
    }
}
