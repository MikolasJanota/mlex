/*
 * File:  BinaryFunction.h
 * Author:  mikolas
 * Created on:  Fri Aug 7 17:05:36 WEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#pragma once
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

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
        for (size_t i = 0; i < _order; i++) {
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

class BinaryFunctionPtr_less {
  public:
    inline bool operator()(const BinaryFunction *s1,
                           const BinaryFunction *s2) const {
        return s1->is_less(*s2);
    }
};

class BinaryFunctionPtr_equal {
  public:
    inline bool operator()(const BinaryFunction *s1,
                           const BinaryFunction *s2) const {
        return s1->is_equal(*s2);
    }
};

class BinaryFunctionPtr_hash {
  public:
    inline size_t operator()(const BinaryFunction *s) const {
        return s->get_hash();
    }
};
