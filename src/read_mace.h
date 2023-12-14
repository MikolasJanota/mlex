/*
 * File:  read_mace.h
 * Author:  mikolas
 * Created on:  Wed Dec 21 13:18:39 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once

#include <zlib.h>

#include <memory>
#include <vector>

#include "binary_function.h"
#include "reader.h"
#include "options.h"

class ReadMace {
  public:
    ReadMace(Output &output, gzFile &input_file);
    size_t read(int max);
    void clear() { d_functions.clear(); }
    const std::vector<std::unique_ptr<BinaryFunction>> &functions() const {
        return d_functions;
    }

  private:
    Output &d_output;
    gzFile &d_input_file;
    Reader d_in;
    std::vector<std::unique_ptr<BinaryFunction>> d_functions;
};
