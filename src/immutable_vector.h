/*
 * File:  ImmutableVector.h
 * Author:  mikolas
 * Created on:  Fri, May 15, 2015 4:47:20 PM
 * Copyright (C) 2015, Mikolas Janota
 */
#pragma once
#include <cassert>
#include <iostream>
#include <iterator>
#include <vector>

template <class T, class HashFun, class Eq> class const_ImmutableVectorIterator;

template <class T, class HashFun = std::hash<T>, class Eq = std::equal_to<T>>
class [[nodiscard]] ImmutableVector {
  public:
    typedef const_ImmutableVectorIterator<T, HashFun, Eq> const_iterator;

  public:
    inline ImmutableVector() {
        _data = nullptr;
        _size = 0;
        _hash_code = EMPTY_HASH;
    }

    explicit ImmutableVector(const std::vector<T> &es);
    virtual ~ImmutableVector() { decrease(); }

    inline ImmutableVector(const ImmutableVector<T, HashFun, Eq> &ls)
        : _hash_code(ls._hash_code), _size(ls._size), _data(ls._data) {
        if (_data != nullptr)
            ++(_data->count);
    }

    ImmutableVector<T, HashFun, Eq> &
    operator=(const ImmutableVector<T, HashFun, Eq> &ls) {
        decrease();
        _hash_code = ls._hash_code;
        _size = ls._size;
        _data = ls._data;
        if (_data != nullptr)
            ++(_data->count);
        return *this;
    }

    bool equals(const ImmutableVector<T, HashFun, Eq> &other) const;
    std::ostream &print(std::ostream &out) const;
    inline size_t hash_code() const { return _hash_code; }
    inline size_t size() const { return _size; }
    inline bool is_empty() const { return _size == 0; }
    inline const_iterator begin() const;
    inline const_iterator end() const;

    inline const T operator[](size_t index) const {
        assert(index < _size);
        return _data->elements[index];
    }

  private:
    struct Data {
        size_t count;
        T *elements;
    };
    static const size_t EMPTY_HASH = 3889;
    size_t _hash_code;
    size_t _size;
    Data *_data;
    inline size_t decrease(); // decrease reference counter
};

template <class T, class HashFun, class Eq>
inline size_t ImmutableVector<T, HashFun, Eq>::decrease() {
    if (_data == nullptr)
        return 0;
    assert(_data->count);
    const size_t nv = --(_data->count);
    if ((_data->count) == 0) {
        delete[] _data->elements;
        delete _data;
    }
    _data = nullptr;
    return nv;
}

template <class T, class HashFun, class Eq>
bool ImmutableVector<T, HashFun, Eq>::equals(
    const ImmutableVector<T, HashFun, Eq> &other) const {
    if (other._size != _size)
        return false;
    if (other._data == _data)
        return true;
    for (size_t i = 0; i < _size; ++i)
        if (!Eq()(_data->elements[i], other._data->elements[i]))
            return false;
    return true;
}

template <class T, class HashFun, class Eq>
class const_ImmutableVectorIterator
    : public std::iterator<std::forward_iterator_tag, T> {
  public:
    const_ImmutableVectorIterator(const ImmutableVector<T, HashFun, Eq> &ls,
                                  size_t x)
        : d_ls(ls), d_i(x) {}
    const_ImmutableVectorIterator(
        const const_ImmutableVectorIterator<T, HashFun, Eq> &mit)
        : d_ls(mit.d_ls), d_i(mit.d_i) {}
    const_ImmutableVectorIterator &operator++() {
        ++d_i;
        return *this;
    }

    const_ImmutableVectorIterator &
    operator=(const const_ImmutableVectorIterator<T, HashFun, Eq> &rhs) {
        assert(&d_ls == &(rhs.d_ls));
        d_i = rhs.d_i;
        return *this;
    }

    bool operator==(const const_ImmutableVectorIterator<T, HashFun, Eq> &rhs) {
        assert(&d_ls == &(rhs.d_ls));
        return d_i == rhs.d_i;
    }
    bool operator!=(const const_ImmutableVectorIterator<T, HashFun, Eq> &rhs) {
        assert(&d_ls == &(rhs.d_ls));
        return d_i != rhs.d_i;
    }
    const T operator*() const { return d_ls[d_i]; }

  private:
    const ImmutableVector<T, HashFun, Eq> &d_ls;
    size_t d_i;
};

template <class T, class HashFun, class Eq>
inline typename ImmutableVector<T, HashFun, Eq>::const_iterator
ImmutableVector<T, HashFun, Eq>::begin() const {
    return const_ImmutableVectorIterator<T, HashFun, Eq>(*this, 0);
}

template <class T, class HashFun, class Eq>
inline typename ImmutableVector<T, HashFun, Eq>::const_iterator
ImmutableVector<T, HashFun, Eq>::end() const {
    return const_ImmutableVectorIterator<T, HashFun, Eq>(*this, _size);
}

template <class T, class HashFun, class Eq>
ImmutableVector<T, HashFun, Eq>::ImmutableVector(const std::vector<T> &es) {
    const size_t vsz = es.size();
    if (vsz == 0) {
        _data = nullptr;
        _size = 0;
        _hash_code = EMPTY_HASH;
        return;
    }
    _data = new Data();
    _data->elements = new T[vsz];
    _data->count = 1;
    T *const elements = _data->elements;
    _size = vsz;
    _hash_code = 7;
    for (size_t i = _size; i--;) {
        elements[i] = es[i];
        _hash_code = _hash_code * 31 + HashFun()(elements[i]);
    }
}

template <class T, class HashFun = std::hash<T>, class Eq = std::equal_to<T>>
class ImmutableVector_equal {
  public:
    inline bool operator()(const ImmutableVector<T, HashFun, Eq> &v1,
                           const ImmutableVector<T, HashFun, Eq> &v2) const {
        return v1.equals(v2);
    }
};

template <class T, class HashFun = std::hash<T>, class Eq = std::equal_to<T>>
struct ImmutableVector_hash {
    inline size_t operator()(const ImmutableVector<T, HashFun, Eq> &ls) const {
        return ls.hash_code();
    }
};

template <class T, class HashFun, class Eq>
std::ostream &operator<<(std::ostream &o,
                         const ImmutableVector<T, HashFun, Eq> &es) {
    o << "[";
    for (size_t i = 0; i < es.size(); ++i)
        o << (i ? " " : "") << es[i];
    return o << "]";
}
