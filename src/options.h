/*
 * File:  options.h
 * Author:  mikolas
 * Created on:  Tue Nov 9 12:11:53 CET 2021
 * Copyright (C) 2021, Mikolas Janota
 */
#pragma once
#include "statistics.h"
#include <cstddef>
#include <iostream>
#include <string>

class null_out_buf : public std::streambuf {
  public:
    virtual std::streamsize xsputn(const char *, std::streamsize n) {
        return n;
    }
    virtual int overflow(int) { return 1; }
};

class null_out_stream : public std::ostream {
  public:
    null_out_stream() : std::ostream(&buf) {}

  private:
    null_out_buf buf;
};

enum SearchType { lin_us, lin_su, bin, bin2, adaptive };
enum ExplicitSearchType { lowering, left_to_right, binary };

struct Options {
    std::string file_name;
    int verbose = 0;
    bool print = false;
    bool graph = false;
    bool explicit_solver = false;
    bool mace_format = false;
    bool unique = false;
    bool color = false;
    bool opt1stRow = false;
    bool simp_sat_row = false;
    bool budgeting = false;
    bool id_elements = false;
    bool invariants = false;
    bool last_solution = false;
    bool use_hash_table = false;
    bool diagonal = false;
    bool budget_idem = false;
    bool inv_ord = false;
    bool distance_invariant = false;
    size_t seq_counter_lits = 10;
    SearchType search_type = SearchType::lin_us;
    ExplicitSearchType explicit_search_type = ExplicitSearchType::lowering;
    const char *comment_prefix = "#";
    std::string diag_file;
};

class Output {
  public:
    null_out_stream cnul; // My null stream.

    Output(Options &options, StatisticsManager &statistics)
        : d_options(options), d_statistics(statistics) {}

    inline std::ostream &ccomment(int level = 0) {
        if (level <= d_options.verbose)
            return std::cout;
        else
            return cnul;
    };
    inline std::ostream &comment(int level = 0) {
        if (d_options.verbose >= level)
            return (std::cout << d_options.comment_prefix << " ").flush();
        else
            return cnul;
    };
    const Options &d_options;
    StatisticsManager &d_statistics;
};
