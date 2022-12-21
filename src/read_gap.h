/*
 * File:  ReadGAP.hh
 * Author:  mikolas
 * Created on:  Fri Oct 2 09:31:16 CEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#pragma once
#include<zlib.h>
#include<vector>
#include<memory>
#include"binary_function.h"
#include"fmtutils.hh"

class ReadGAP {
    public:
        explicit ReadGAP(gzFile& input_file);
        virtual ~ReadGAP() = default;
        void read();
        bool has_f() const { return static_cast<bool>(_f); }
        const BinaryFunction& f() { assert(has_f()); return *_f; };
    private:
        gzFile _input_file;
        std::unique_ptr<BinaryFunction> _f;
};
