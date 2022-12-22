/*
 * File:  trie.h
 * Author:  mikolas
 * Created on:  Thu Dec 22 16:58:34 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once

#include "binary_function.h"
#include <algorithm>
#include <utility>
#include <vector>

/** A single node of the IndexTrie. */
struct ModelTrieNode {
    std::vector<ModelTrieNode *> d_children;
    ModelTrieNode *              d_blank = nullptr;
};

/** */
class ModelTrie {
  public:
    /*  Construct the trie */
    ModelTrie() : d_root{new ModelTrieNode()} {}

    virtual ~ModelTrie() { freeRec(d_root); }

    bool add(const BinaryFunction &m);

  private:
    /**  the root of the trie, becomes null, if all tuples should match */
    ModelTrieNode *d_root;

    /** Auxiliary recursive function for cleanup. */
    void freeRec(ModelTrieNode *n);

    /** Add master values  starting from index  to a given subtree. The
     * cardinality represents the number of non-blank elements left. */
    std::pair<bool, ModelTrieNode *> addRec(ModelTrieNode *n, size_t index,
                                            const BinaryFunction &m);
};
