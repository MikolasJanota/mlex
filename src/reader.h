/*
 * File:   Reader.hh
 * Author: mikolas
 *
 * Created on January 12, 2011, 4:19 PM
 */
#pragma once
#include "stream_buffer.h"
#include <cstdio>
#include <iostream>
#include <zlib.h>
class Reader {
  public:
    Reader(gzFile &zf);
    virtual ~Reader();
    [[nodiscard]] int operator*() const { return cur; }
    void operator++() { update_cur(); }
    void skip_whitespace();
    inline size_t get_line_number() const { return lnn; }

    void skip_line();
    void update_cur();

    int parse_int() {
        skip_whitespace();
        auto &in = *this;
        const bool neg = *in == '-';
        if (*in == '-' || *in == '+')
            ++in;
        if (*in < '0' || *in > '9') {
            std::cerr << get_line_number() << ":expected digit instead of '"
                      << static_cast<char>(*in) << "'" << std::endl;
            exit(EXIT_FAILURE);
        }
        int val = 0;
        for (; *in >= '0' && *in <= '9'; ++in)
            val = val * 10 + (*in - '0');
        return neg ? -val : val;
    }

  private:
    size_t lnn;
    StreamBuffer buf;
    int cur;
};

