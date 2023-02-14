/*
 * File:  ReadGAP.hh
 * Author:  mikolas
 * Created on:  Fri Oct 2 09:31:16 CEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#pragma once
#include "binary_function.h"
#include "fmtutils.hh"
#include "options.h"
#include <memory>
#include <vector>
#include <zlib.h>

class ReadGAP {
  public:
    explicit ReadGAP(Output &output, gzFile &input_file);
    virtual ~ReadGAP() = default;
    size_t read(int max);
    void clear() { d_functions.clear(); }

    const std::vector<std::unique_ptr<BinaryFunction>> &functions() const {
        return d_functions;
    }

  private:
    Output &d_output;
    gzFile &d_input_file;
    StreamBuffer d_in;
    bool d_open = false;
    bool d_closed = false;
    std::vector<std::unique_ptr<BinaryFunction>> d_functions;
    std::unique_ptr<BinaryFunction> read_fun();
};
