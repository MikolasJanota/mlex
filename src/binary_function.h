/*
 * File:  BinaryFunction.h
 * Author:  mikolas
 * Created on:  Fri Aug 7 17:05:36 WEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#pragma once
#include "CLI11.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <math.h>
#include <vector>

inline size_t calculate_bits_needed(size_t order) {
    assert(order > 0);
    size_t rv = 0;
    for (--order; order > 0; order >>= 1)
        rv++;
    return rv;
}

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

    void print_mace(std::ostream &output,
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
            size_t val = d_f.get(d_pos++);
            /* std::cerr<<"v:"<<val<<'\n'; */
            for (auto i = d_fits; i--;) {
                d_buf.push_back(val & d_mask);
                /* std::cerr<<" p:"<<d_buf.back()<<" "; */
                val >>= d_bits;
            }
            std::reverse(d_buf.begin(), d_buf.end());
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

    CompFunction make() {
        if (!d_buf.empty())
            flush();
        const auto rv = CompFunction(d_order, d_arity, d_data_size, d_data);
        d_data = nullptr;
        return rv;
    }

    void reset() {
        assert(d_buf.empty());
        d_pos = 0;
        d_buf.clear();
        if (d_data == nullptr)
            d_data = new uint64_t[d_data_size];
    }

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

class BinaryFunction {
  public:
    explicit BinaryFunction(size_t order)
        : _order(order), _values(order * order, -1) {}

    void set(size_t i, size_t j, size_t val) {
        assert(i < _order && j < _order);
        _values[i * _order + j] = val;
    }

    size_t get(size_t i, size_t j) const {
        assert(i < _order && j < _order);
        return _values[i * _order + j];
    }

    size_t order() const { return _order; }

    bool is_equal(const BinaryFunction &other) const {
        return _order == other._order && _values == other._values;
    }

    bool is_less(const BinaryFunction &other) const {
        if (_order != other._order)
            return _order < other._order;
        const auto sz = _order * _order;
        for (size_t i = 0; i < sz; i++) {
            if (_values[i] != other._values[i])
                return _values[i] < other._values[i];
        }
        return false;
    }

    inline uint64_t get_hash() const {
        assert(_hash_defined);
        return _hash;
    }

    void setup_hash() {
        _hash = 7 * _order;
        for (const auto &val : _values) {
            uint64_t x = static_cast<uint64_t>(val);
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            _hash ^= x + 0x9e3779b9 + (_hash << 6) + (_hash >> 2);
        }
        _hash_defined = true;
    }

    void print(std::ostream &out) const;
    void set_name(const std::string &s) { _name = s; }
    const std::string &get_name() const { return _name; }
    void set_additional_info(const std::string &s) { _additional_info = s; }
    const std::string &get_additional_info() const { return _additional_info; }

    void print_gap(std::ostream &output);
    void print_mace(std::ostream &output);

  private:
    const size_t _order;
    bool _hash_defined = false;
    uint64_t _hash;
    std::string _name;
    std::string _additional_info;
    std::vector<size_t> _values;
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
