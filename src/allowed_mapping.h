/*
 * File:  allowed_mapping.h
 * Author:  mikolas
 * Created on:  Fri Aug 11 15:48:46 CEST 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include <cstddef>
#include <vector>

class AllowedMapping {
  public:
    AllowedMapping(size_t n)
        : d_n(n), d_allowed(d_n, std::vector<bool>(d_n, true)) {}

    bool allowed(size_t from, size_t to) const { return d_allowed[from][to]; }

    const std::vector<bool> &get_allowed(size_t from) const {
        return d_allowed[from];
    }

    bool disallow(size_t from, size_t to) {
        const auto oldv = d_allowed[from][to];
        d_allowed[from][to] = false;
        return oldv;
    }

  protected:
    const size_t d_n;
    std::vector<std::vector<bool>> d_allowed;
};
