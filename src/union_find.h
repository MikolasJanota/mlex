/*
 * File:  union_find.h
 * Author:  mikolas
 * Created on:  13 Jan 2018 17:38:22
 * Copyright (C) 2018, Mikolas Janota
 */
#pragma once
#include "auxiliary.h"
#include <unordered_map>
template <class T, class HashFun = std::hash<T>, class Eq = std::equal_to<T>>
class UnionFind {
  public:
    bool find(const T &a, T &repr) const {
        const auto i = _parent.find(a);
        if (i == _parent.end())
            return false;
        repr = get_repr(a, i->second);
        return true;
    }

    T get_repr(const T &a) const {
        T repr;
        VERIFY(find(a, repr));
        return repr;
    }

    T find_or_make(const T &a) {
        const auto i = _parent.find(a);
        if (i == _parent.end()) {
            _parent.insert(i, {a, a});
            _rank[a] = 0;
            return a;
        } else {
            return get_repr(a, i->second);
        }
    }

    void unify(const T &a, const T &b) {
        const T ar = find_or_make(a);
        const T br = find_or_make(b);
        if (ar == br)
            return;
        const auto arn = rank(ar);
        const auto brn = rank(br);
        if (arn < brn)
            _parent[ar] = br;
        else if (arn > brn)
            _parent[br] = ar;
        else {
            _parent[ar] = br;
            _rank[br] = brn + 1;
        }
    }

    bool same(const T &a, const T &b) const {
        T ra, rb;
        return find(a, ra) && find(b, rb) && ra == rb;
    }

  private:
    std::unordered_map<T, T, HashFun, Eq> _parent;
    std::unordered_map<T, int, HashFun, Eq> _rank;

    int rank(const T &a) const {
        const auto i = _rank.find(a);
        assert(i != _rank.end());
        return i->second;
    }

    T get_parent(const T &a) const {
        const auto i = _parent.find(a);
        assert(i != _parent.end());
        return i->second;
    }

    T get_repr(T a, T parent) const {
        while (a != parent) {
            a = parent;
            parent = get_parent(a);
        }
        return a;
    }
};
