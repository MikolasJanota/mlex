/*
 * File:  read_diags.h
 * Author:  mikolas
 * Created on:  Mon Feb 13 08:38:43 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once

#include <zlib.h>

#include <memory>
#include <vector>

#include "binary_function.h"
#include "fmtutils.hh"
#include "options.h"

class ReadDiags {
  public:
    ReadDiags(Output &output, gzFile &input_file);
    size_t read(int max);
    void clear() { d_diags.clear(); }
    const std::vector<std::vector<size_t>> &diags() const { return d_diags; }

  private:
    Output &d_output;
    gzFile &d_input_file;
    StreamBuffer d_in;
    std::vector<std::vector<size_t>> d_diags;
    bool d_open = false;
    bool d_closed = false;
    size_t d_read_diags = 0;
};
