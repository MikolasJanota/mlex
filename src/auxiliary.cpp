/*
 * File:  auxiliary.cc
 * Author:  mikolas
 * Created on:  Wed, Dec 03, 2014 4:00:16 PM
 * Copyright (C) 2014, Mikolas Janota
 */
#include <cassert>
#include <iostream>
#include <vector>
std::ostream &show_permutation(std::ostream &out, std::vector<size_t> &perm) {
    const auto n = perm.size();
    std::vector<bool> visited(n, false);
    while (1) {
        size_t s = 0;
        while (visited[s] && s < n)
            s++;
        if (s == n)
            return out;
        out << "(" << s << std::flush;
        visited[s] = true;
        for (auto i = perm[s]; i != s; i = perm[i]) {
            assert(i < n);
            out << " " << i << std::flush;
            visited[i] = true;
        }
        out << ")";
    }
}

