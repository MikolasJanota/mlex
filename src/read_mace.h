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
#include "fmtutils.hh"

class ReadMace {
  public:
    explicit ReadMace(gzFile &input_file);
    void                                                read();
    const std::vector<std::unique_ptr<BinaryFunction>> &functions() const {
        return d_functions;
    }

  private:
    gzFile                                       d_input_file;
    std::vector<std::unique_ptr<BinaryFunction>> d_functions;
};
