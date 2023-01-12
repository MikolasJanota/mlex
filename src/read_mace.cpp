/*
 * File:  read_mace.cpp
 * Author:  mikolas
 * Created on:  Wed Dec 21 13:18:43 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */

#include "read_mace.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include "fmtutils.hh"

ReadMace::ReadMace(Output &output, gzFile &input_file)
    : d_output(output), d_input_file(input_file), d_in(d_input_file) {}

static void skip(StreamBuffer &sb) {
    while (1) {
        skipWhitespace(sb);
        if (*sb == '%')
            skipLine(sb);
        else
            break;
    }
}

static void match_string(StreamBuffer &sb, const char *s) {
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

static void match_char(StreamBuffer &sb, char c) {
    skip(sb);
    const char rc = *sb;
    if (rc != c) {
        std::cerr << c << " expected instead of '" << rc << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
    ++sb;
}

static void match_chars(StreamBuffer &sb, const char *s) {
    while (*s != '\0') {
        match_char(sb, *s);
        s++;
    }
}

size_t ReadMace::read(int max) {
    assert(max >= 0);

    size_t counter = 0;
    while (skip(d_in), *d_in != EOF && max--) {
        match_string(d_in, "interpretation");
        match_char(d_in, '(');
        const auto i = parseInt(d_in);
        if (i < 1) {
            std::cerr << "invalid order " << i << std::endl;
            exit(EXIT_FAILURE);
        }
        const size_t order = static_cast<size_t>(i);
        d_functions.push_back(std::make_unique<BinaryFunction>(order));
        auto &f = d_functions.back();

        match_char(d_in, ',');
        match_char(d_in, '[');
        std::stringstream buf1;
        while (*d_in != EOF && *d_in != ']') {
            buf1 << static_cast<char>(*d_in);
            ++d_in;
        }
        f->set_additional_info(buf1.str());
        match_char(d_in, ']');
        match_char(d_in, ',');
        match_char(d_in, '[');
        match_string(d_in, "function");
        match_char(d_in, '(');
        std::stringstream buf2;
        while (*d_in != EOF && *d_in != '(') {
            buf2 << *d_in;
            ++d_in;
        }
        f->set_name(buf2.str());
        match_chars(d_in, "(_,_),[");
        for (size_t i = 0; i < order; i++) {
            for (size_t j = 0; j < order; j++) {
                f->set(i, j, parseInt(d_in));
                if (i + 1 < order || j + 1 < order)
                    match_char(d_in, ',');
            }
        }
        match_chars(d_in, "])]).");
        counter++;
        d_output.comment(1)
            << "read function " << f->get_name() << " order:" << f->order()
            << " " << f->get_additional_info() << std::endl;
    }
    return counter;
}
