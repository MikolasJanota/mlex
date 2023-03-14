/*
 * File:  fixed_elem_map.h
 * Author:  mikolas
 * Created on:  Tue Mar 14 16:42:25 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once

#include "options.h"
#include <cassert>
#include <cstddef>
#include <set>

class FixedElemMap {
  public:
    FixedElemMap(Output &output, size_t order)
        : d_output(output), d_order(order), d_src2dst(d_order, 2 * order + 1),
          d_dst2src(d_order, 2 * order + 1) {}
    virtual ~FixedElemMap() {}

    size_t fixed_count() const { return d_fixed_src_elements.size(); }

    const std::vector<size_t> &fixed_src_elements() const {
        return d_fixed_src_elements;
    }

    bool is_fixed_src(size_t src) const { return d_src2dst[src] < d_order; }
    bool is_fixed_dst(size_t dst) const { return d_dst2src[dst] < d_order; }

    size_t src2dst(size_t src) const {
        assert(is_fixed_src(src));
        return d_src2dst[src];
    }

    size_t dst2src(size_t dst) const {
        assert(is_fixed_dst(dst));
        return d_dst2src[dst];
    }

    bool set(size_t src, size_t dst) {
        assert(!is_fixed_src(src) || src2dst(src) == dst);
        assert(!is_fixed_src(src) || dst2src(dst) == src);
        if (is_fixed_src(src))
            return false;
        assert(!is_fixed_dst(dst));
        d_src2dst[src] = dst;
        d_dst2src[dst] = src;
        d_output.d_statistics.fixedElements->inc();
        d_fixed_src_elements.push_back(src);
        d_fixed_count++;
        return true;
    }

  private:
    Output &d_output;
    const size_t d_order;
    size_t d_fixed_count;
    std::vector<size_t> d_src2dst;
    std::vector<size_t> d_dst2src;
    std::vector<size_t> d_fixed_src_elements;
};
