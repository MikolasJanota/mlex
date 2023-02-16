/*
 * File:   auxiliary.hh
 * Author: mikolas
 *
 * Created on October 12, 2011
 */
#pragma once
#include <set>
#include <sys/time.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#ifndef __MINGW32__
#include <sys/resource.h>
#endif
#include <cassert>
#include <iomanip>
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

#define SHOW_TIME1(t) std::fixed << std::setprecision(1) << (t)
#define SHOW_TIME2(t) std::fixed << std::setprecision(2) << (t)
#define SHOW_TIME4(t) std::fixed << std::setprecision(4) << (t)
#define SHOW_TIME0(t) std::fixed << std::setprecision(0) << (t)

#ifdef NDEBUG
#define SHOW_TIME(t) SHOW_TIME4(t)
#else
#define SHOW_TIME(t) SHOW_TIME1(t)
#endif

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

template <class K>
inline bool insert_chk(std::unordered_set<K> &es, const K &e) {
    const auto i = es.insert(e);
    const bool a = i.second;
    assert(a);
    return a;
}

template <class M, class K> bool contains(const M &m, const K &e) {
    return m.find(e) != m.end();
}

template <class K> const K &first(const std::set<K> &m) {
    assert(!m.empty());
    return *(m.begin());
}

template <class C> std::ostream &print_set(std::ostream &out, const C &m) {
    out << "{";
    for (const auto &e : m)
        out << " " << e;
    return out << " }";
}
template <class C>
std::ostream &print_vec(std::ostream &out, const C &m, size_t offset = 0) {
    out << "[";
    bool f = true;
    for (const auto &e : m) {
        out << (f ? "" : ",") << (offset + e);
        f = false;
    }
    return out << "]";
}
