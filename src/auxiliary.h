/*
 * File:   auxiliary.hh
 * Author: mikolas
 *
 * Created on October 12, 2011
 */
#pragma once
#include <sys/time.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <assert.h>
#include <iomanip>
#include <string.h>
#define OUT std::cout

#define __PL (std::cerr << __FILE__ << ":" << __LINE__ << std::endl).flush();

#define SATSPC Minisat

#define VERIFY(condition)                                                      \
    do {                                                                       \
        if (condition)                                                         \
            ;                                                                  \
        else                                                                   \
            assert(0);                                                         \
    } while (0)

#ifdef __MINGW32__
inline double read_cpu_time() { return 0; }
#else
inline double read_cpu_time() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}
#endif

#define SHOW_TIME(t)  std::fixed << std::setprecision(4) << (t)
#define SHOW_TIME2(t) std::fixed << std::setprecision(2) << (t)
#define SHOW_TIME0(t) std::fixed << std::setprecision(0) << (t)

template <class K, class V>
bool contains(const std::unordered_map<K, V> &es, const K &e) {
    return es.find(e) != es.end();
}

template <class K, class V>
V get(const std::unordered_map<K, V> &es, const K &e) {
    const auto j = es.find(e);
    assert(j != es.end());
    return j->second;
}

inline bool contains(const std::vector<bool> &es, int e) {
    assert(e >= 0);
    const auto ix = (size_t)e;
    return ix < es.size() ? es[ix] : false;
}

template <class K>
inline bool insert_chk(std::unordered_set<K> &es, const K &e) {
    const auto i = es.insert(e);
    const bool a = i.second;
    assert(a);
    return a;
}

inline bool erase(std::vector<bool> &es, int e) {
    assert(e >= 0);
    const auto ix = (size_t)e;
    if (ix >= es.size())
        return false;
    const bool rv = es[ix];
    es[ix]        = false;
    return rv;
}

inline bool insert(std::vector<bool> &es, int e) {
    assert(e >= 0);
    const auto ix = (size_t)e;
    if (ix >= es.size())
        es.resize(ix + 1, false);
    const bool rv = es[ix];
    es[ix]        = e;
    return rv;
}

template <class K> bool contains(const std::unordered_set<K> &es, const K &e) {
    return es.find(e) != es.end();
}
