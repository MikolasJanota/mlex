/*
 * File:  BinaryFunction.h
 * Author:  mikolas
 * Created on:  Fri Aug 7 17:05:36 WEST 2020
 * Copyright (C) 2020, Mikolas Janota
 */
#pragma once
#include <cassert>
#include <iostream>
#include <cstddef> // for size_t
#include <string>   // for string
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

    void print(std::ostream &out) const;
    void set_name(const std::string &s) { _name = s; }
    const std::string &get_name() const { return _name; }
    void set_additional_info(const std::string &s) { _additional_info = s; }
    const std::string &get_additional_info() const { return _additional_info; }

    void print_gap(std::ostream &output);
    void print_mace(std::ostream &output);

  private:
    const size_t _order;
    std::string _name;
    std::string _additional_info;
    std::vector<size_t> _values;
};
