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

ReadMace::ReadMace(Output &output, gzFile &input_file)
    : d_output(output), d_input_file(input_file), d_in(d_input_file) {}

static void skip(Reader &sb) {
    while (1) {
        sb.skip_whitespace();
        if (*sb == '%')
            sb.skip_line();
        else
            break;
    }
}

static int parse_int(Reader &sb) {
    skip(sb);
    return sb.parse_int();
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

static void match_char(Reader &sb, char c) {
    skip(sb);
    const char rc = *sb;
    if (rc != c) {
        std::cerr << sb.get_line_number() << ":expected '" << c
                  << "' instead of '" << static_cast<char>(*sb) << "'"
                  << std::endl;
        std::cerr << "around: \"";
        for (size_t i = 0; i < 10 && *sb != EOF; i++, ++sb)
            std::cerr << static_cast<char>(*sb);
        std::cerr << "\"" << std::endl;
        exit(EXIT_FAILURE);
    }
    ++sb;
}

static void match_chars(Reader &sb, const char *s) {
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
        const auto i = parse_int(d_in);
        if (i < 1) {
            std::cerr << "invalid order " << i << std::endl;
            exit(EXIT_FAILURE);
        }
        const size_t order = static_cast<size_t>(i);
        d_functions.push_back(std::make_unique<BinaryFunction>(order));
        auto &f = d_functions.back();

        match_char(d_in, ',');
        match_char(d_in, '[');
        {
            std::stringstream buf1;
            while (*d_in != EOF && *d_in != ']') {
                buf1 << static_cast<char>(*d_in);
                ++d_in;
            }
            f->set_additional_info(buf1.str());
        }
        match_char(d_in, ']');
        match_char(d_in, ',');
        match_char(d_in, '[');
        match_string(d_in, "function");
        match_char(d_in, '(');
        {
            std::stringstream buf2;
            while (*d_in != EOF && *d_in != '(') {
                buf2 << static_cast<char>(*d_in);
                ++d_in;
            }
            f->set_name(buf2.str());
        }
        match_chars(d_in, "(_,_),[");
        for (size_t i = 0; i < order; i++) {
            for (size_t j = 0; j < order; j++) {
                const auto val = parse_int(d_in);
                if (val < -1 || val >= static_cast<int>(order)) {
                    std::cerr << "value '" << val << "' out of range"
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (val >= 0)
                  f->set(i, j, val);
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
