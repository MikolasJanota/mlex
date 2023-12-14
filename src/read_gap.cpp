/*
 * File:  ReadGAP.cpp
 * Author:  mikolas
 * Created on:  Fri Oct 2 09:31:24 CEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#include "read_gap.h"
#include "auxiliary.h"
#include <cstdlib>
#include <iostream>
#include <vector>
using std::vector;

ReadGAP::ReadGAP(Output &output, gzFile &input_file)
    : d_output(output), d_input_file(input_file), d_in(d_input_file) {}

static void skip(Reader &sb) {
    for (;;) {
        sb.skip_whitespace();
        if (*sb == '#')
            sb.skip_line();
        else
            break;
    }
}

static void match_char(Reader &sb, char c) {
    skip(sb);
    if (*sb != c) {
        std::cerr << sb.get_line_number() << ":expected '" << c
                  << "' instead of '" << static_cast<char>(*sb) << "'"
                  << std::endl;
        std::cerr << "around: \"";
        const auto fc = *sb;
        for (size_t i = 0; i < 10 && *sb != EOF; i++, ++sb)
            std::cerr << static_cast<char>(*sb);
        std::cerr << "\"" << std::endl;
        if (fc == 'i')
            std::cerr << "NOTE: for MACE format, the -m option is needed."
                      << std::endl;
        exit(EXIT_FAILURE);
    }
    ++sb;
}

static int check_next_char(Reader &sb) {
    skip(sb);
    return *sb;
}

static void match_string(Reader &sb, const char *s) {
    skip(sb);
    const auto olds = s;
    for (; *s; ++sb, s++) {
        if (*sb == EOF) {
            std::cerr << sb.get_line_number()
                      << ":End of file when looking for '" << olds << "'."
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        const char rc = *sb;
        if (rc != *s) {
            std::cerr << sb.get_line_number() << ":Unexpected character '" << rc
                      << "' when looking for '" << *s << "' in '" << olds
                      << "'." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

size_t ReadGAP::read(int max) {
    if (max <= 0 || d_closed)
        return 0;
    if (!d_open) {
        if (check_next_char(d_in) == 't') {
            match_string(d_in, "tables");
            match_string(d_in, ":=");
        }
        match_char(d_in, '[');
        d_open = true;
    }
    int cnt = 0;
    while (cnt < max && check_next_char(d_in) == '[') {
        d_functions.push_back(read_fun());
        cnt++;
    }
    if (check_next_char(d_in) == ']') {
        match_char(d_in, ']');
        d_closed = true;
    } else
        match_char(d_in, ',');
    return cnt;
}

std::unique_ptr<BinaryFunction> ReadGAP::read_fun() {
    match_char(d_in, '['); // open table

    vector<int> row;
    vector<vector<int>> rows;
    size_t order = -1;

    while (true) {
        match_char(d_in, '['); // open row
        row.clear();
        while (true) {
            skip(d_in);
            row.push_back(d_in.parse_int());
            if (check_next_char(d_in) == ']') {
                ++d_in; // close row
                break;
            }
            match_char(d_in, ',');
        }
        if (rows.size() && row.size() != order) {
            std::cerr << "invalid row size " << std::endl;
            exit(EXIT_FAILURE);
        }
        order = row.size();
        rows.push_back(row);
        if (check_next_char(d_in) == ']') {
            ++d_in; // close tables
            break;
        }
        match_char(d_in, ',');
    }
    if (rows.size() != order) {
        std::cerr << "invalid number of rows" << std::endl;
        exit(EXIT_FAILURE);
    }
    auto rv = std::make_unique<BinaryFunction>(order);
    for (size_t i = 0; i < order; i++) {
        for (size_t j = 0; j < order; j++) {
            const auto number = rows[i][j];
            if (number <= 0 || number > static_cast<int>(order)) {
                std::cerr << "invalid number " << number << std::endl;
                exit(EXIT_FAILURE);
            }
            rv->set(i, j, number - 1);
        }
    }
    return rv;
}
