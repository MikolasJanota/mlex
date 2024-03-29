/*
 * File:  comp_function.h
 * Author:  mikolas
 * Created on:  Fri Jan 27 15:39:42 CET 2023
 * Copyright (C) 2023, Mikolas Janota
 */
#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

inline size_t calculate_bits_needed(size_t order) {
    assert(order > 0);
    size_t rv = 0;
    for (--order; order > 0; order >>= 1)
        rv++;
    return rv;
}

/* Calculate how many values fit in 64-bits. */
inline size_t calculate_fits(size_t order) {
    assert(order > 1);
    return std::numeric_limits<uint64_t>::digits / calculate_bits_needed(order);
}

inline size_t calculate_cell_count(size_t order, size_t arity) {
    assert(order > 1);
    size_t sz = 1;
    while (arity--)
        sz *= order;
    return sz;
}

inline size_t calculate_data_size(size_t order, size_t arity) {
    assert(order > 1);
    const auto sz = calculate_cell_count(order, arity);
    const auto f = calculate_fits(order);
    return 1 + (sz - 1) / f;
}

/*A representation of a multiplication table that's meant to occupy little
 * space. Multiple values are masked into 64bits.*/
class CompFunction {
  public:
    CompFunction(size_t order, size_t arity, size_t data_size, uint64_t *data)
        : d_order(order), d_arity(arity), d_data_size(data_size), d_data(data) {
        assert(data_size == calculate_data_size(d_order, d_arity));
        setup_hash();
    }

    void free() {
        if (d_data != nullptr)
            delete[] d_data;
        d_data = nullptr;
    }

    size_t order() const { return d_order; }
    size_t arity() const { return d_arity; }
    uint64_t get(size_t i) const { return d_data[i]; }
    inline uint64_t get_hash() const { return d_hash; }

    std::ostream &print_gap(std::ostream &output) const;
    std::ostream &print_mace(std::ostream &output,
                             const std::string &additional_info) const;
    void set_name(const std::string &s) { _name = s; }
    const std::string &get_name() const { return _name; }
    uint64_t *get_data() const { return d_data; }

    bool is_equal(const CompFunction &other) const {
        if ((d_order != other.d_order) || (d_arity != other.d_arity))
            return false;
        for (auto i = d_data_size; i--;)
            if (d_data[i] != other.d_data[i])
                return false;
        return true;
    }

    bool is_less(const CompFunction &other) const {
        if (d_order != other.d_order)
            return d_order < other.d_order;
        if (d_arity != other.d_arity)
            return d_arity < other.d_arity;
        assert(d_data_size == other.d_data_size);
        for (size_t i = 0; i < d_data_size; i++) {
            if (d_data[i] != other.d_data[i])
                return d_data[i] < other.d_data[i];
        }
        return false;
    }

  private:
    const size_t d_order, d_arity, d_data_size;
    uint64_t *d_data;
    uint64_t d_hash;
    std::string _name;

    void setup_hash() {
        d_hash = 7 * d_order;
        for (size_t i = d_data_size; i--;) {
            uint64_t x = d_data[i];
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            d_hash ^= x + 0x9e3779b9 + (d_hash << 6) + (d_hash >> 2);
        }
    }
};

/* Reads values sequentially from compact representation. */
class CompFunctionReader {
  public:
    CompFunctionReader(const CompFunction &f)
        : d_f(f), d_order(f.order()), d_arity(f.arity()) {
        d_bits = calculate_bits_needed(d_order);
        d_fits = calculate_fits(d_order);
        d_mask = (static_cast<uint64_t>(1) << d_bits) - 1;
    }

    size_t next() {
        if (d_buf.empty()) {
            d_buf.resize(d_fits);
            size_t val = d_f.get(d_pos++);
            for (auto i = d_fits; i--;) {
                d_buf[i] = val & d_mask;
                val >>= d_bits;
            }
        }
        const auto rv = d_buf.back();
        d_buf.pop_back();
        return rv;
    }

  private:
    const CompFunction &d_f;
    const size_t d_order, d_arity;
    std::vector<size_t> d_buf;
    size_t d_pos = 0;
    uint64_t d_mask;
    size_t d_bits = -1, d_fits = -1;
};

/* Helper to build a compact representation. */
class CompFunctionBuilder {
  public:
    CompFunctionBuilder(size_t order, size_t arity)
        : d_order(order), d_arity(arity),
          d_data_size(calculate_data_size(order, arity)) {
        d_bits = calculate_bits_needed(order);
        d_fits = calculate_fits(order);
        reset();
    }

    ~CompFunctionBuilder() {
        if (d_data != nullptr)
            delete[] d_data;
    }

    /* Create a compact function based on whatever has been pushed so far. */
    CompFunction make() {
        if (!d_buf.empty())
            flush();
        const auto rv = CompFunction(d_order, d_arity, d_data_size, d_data);
        d_data = nullptr;
        return rv;
    }

    /* Prepare for a new function. */
    void reset() {
        assert(d_buf.empty());
        d_pos = 0;
        d_buf.clear();
        if (d_data == nullptr)
            d_data = new uint64_t[d_data_size];
    }

    /* Push value into the table. */
    void push(size_t v) {
        assert(v < d_order);
        d_buf.push_back(v);
        assert(d_buf.size() <= d_fits);
        if (d_buf.size() == d_fits)
            flush();
    }

  private:
    const size_t d_order, d_arity, d_data_size;
    std::vector<size_t> d_buf;
    size_t d_pos;
    uint64_t *d_data = nullptr;
    size_t d_bits = -1, d_fits = -1;

    void flush() {
        assert(!d_buf.empty());
        auto &w = d_data[d_pos++];
        w = 0;
        while (!d_buf.empty()) {
            w <<= d_bits;
            w |= d_buf.back();
            d_buf.pop_back();
        }
    }
};

class CompFunction_less {
  public:
    inline bool operator()(const CompFunction &s1,
                           const CompFunction &s2) const {
        return s1.is_less(s2);
    }
};

class CompFunction_equal {
  public:
    inline bool operator()(const CompFunction &s1,
                           const CompFunction &s2) const {
        return s1.is_equal(s2);
    }
};

class CompFunction_hash {
  public:
    inline size_t operator()(const CompFunction &s) const {
        return s.get_hash();
    }
};
