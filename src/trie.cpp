/*
 * File:  trie.cpp
 * Author:  mikolas
 * Created on:  Thu Dec 22 17:03:36 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "trie.h"

void ModelTrie::freeRec(ModelTrieNode *n) {
    if (!n)
        return;
    for (auto c : n->d_children)
        freeRec(c);
    delete n;
}

bool ModelTrie::add(const BinaryFunction &m) {
    const auto     N     = m.order();
    ModelTrieNode *n     = d_root;
    bool           added = false;
    for (size_t index = 0; true; index++) {
        const auto row = index / N;
        if (row >= N)
            return added;
        const auto col = index % N;
        const auto val = m.get(row, col);
        if (n->d_children.empty())
            n->d_children.resize(N, nullptr);

        assert(n->d_children.size() == N);

        if (n->d_children[val] == nullptr) {
            n->d_children[val] = new ModelTrieNode();
            added              = true;
        }
        n = n->d_children[val];
    }
    return added;
}
