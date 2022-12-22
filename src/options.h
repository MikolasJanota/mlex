/*
 * File:  options.h
 * Author:  mikolas
 * Created on:  Tue Nov 9 12:11:53 CET 2021
 * Copyright (C) 2021, Mikolas Janota
 */
#pragma once
#include <cstddef>
#include <string>
struct Options {
    int verbose     = 0;
    int incremental = 0;
    int mace_format = 0;
    int unique      = 0;
};
