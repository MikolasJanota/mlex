/*
 * File:  read_diags.cpp
 * Author:  mikolas
 * Created on:  Mon Feb 13 08:40:27 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */

#include "read_diags.h"
#include "auxiliary.h"

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

ReadDiags::ReadDiags(Output &output, gzFile &input_file)
    : d_output(output), d_input_file(input_file), d_in(d_input_file) {}

static void skip(Reader &sb) { sb.skip_whitespace(); }

static void match_string(Reader &sb, const char *s) {
    skip(sb);
    const auto olds = s;
    for (; *s; ++sb, s++) {
        if (*sb == EOF) {
            std::cerr << "End of file when looking for '" << olds << "'."
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        const char rc = *sb;
        if (rc != *s) {
            std::cerr << "Unexpected character '" << rc
                      << "' when looking for '" << *s << "' in " << olds << "'."
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

static int check_next_char(Reader &sb) {
    skip(sb);
    return *sb;
}

static void match_char(Reader &sb, char c) {
    skip(sb);
    const char rc = *sb;
    if (rc != c) {
        std::cerr << "'" << c << "' expected instead of '" << rc << "'"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    ++sb;
}

size_t ReadDiags::read(int max) {
    assert(max >= 0);
    if (d_closed)
        return 0;
    if (!d_open) {
        if (check_next_char(d_in) == 'm') {
            match_string(d_in, "mindiags");
            match_string(d_in, ":=");
        }
        match_char(d_in, '[');
        d_open = true;
    }

    size_t counter = 0;
    while (skip(d_in), *d_in != EOF && max--) {
        const auto nc = check_next_char(d_in);
        if (nc == ']') {
            match_char(d_in, ']');
            match_char(d_in, ';');
            d_closed = false;
            return counter;
        }
        if (d_read_diags)
            match_char(d_in, ',');
        match_char(d_in, '[');
        // diag
        match_char(d_in, '[');
        std::vector<size_t> diag;
        while (skip(d_in), *d_in != ']') {
            if (!diag.empty())
                match_char(d_in, ',');
            const auto i = d_in.parse_int();
            diag.push_back(i - 1);
        }
        d_diags.push_back(diag);
        match_char(d_in, ']'); // close diag
        match_char(d_in, ',');
        // perm
        while (check_next_char(d_in) == '(') {
            match_char(d_in, '(');
            size_t cyc_size = 0;
            while (check_next_char(d_in) != ')') {
                if (cyc_size++ > 0)
                    match_char(d_in, ',');
                d_in.parse_int();
            }
            match_char(d_in, ')');
        }
        match_char(d_in, ']'); // close diag+perm
        counter++;
        d_read_diags++;
    }
    return counter;
}
