/*
 * File:  trie.cpp
 * Author:  mikolas
 * Created on:  Thu Dec 22 17:03:36 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "trie.h"

bool ModelTrie::add(const BinaryFunction &m) {
    auto [added, root] = addRec(d_root, 0, m);
    d_root             = root;
    return added;
}

void ModelTrie::freeRec(ModelTrieNode *n) {
    if (!n)
        return;
    for (auto c : n->d_children)
        freeRec(c);
    delete n;
}

std::pair<bool, ModelTrieNode *>
ModelTrie::addRec(ModelTrieNode *n, size_t index, const BinaryFunction &m) {
    const auto N   = m.order();
    const auto row = index / N;
    if (row >= N)
        return {false, n};
    const auto col = index % N;
    const auto val = m.get(row, col);

    if (n->d_children.empty())
        n->d_children.resize(N, nullptr);

    assert(n->d_children.size() == N);

    auto old_child = n->d_children[val];
    if (old_child == nullptr) {
        n->d_children[val] = addRec(new ModelTrieNode(), index + 1, m).second;
        return {true, n};
    } else
        return {addRec(old_child, index + 1, m).first, n};
}
