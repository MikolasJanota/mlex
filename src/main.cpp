/*
 * File:  main.cpp
 * Author:  mikolas
 * Created on:  Tue Dec 13 12:04:23 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "CLI11.hpp"
#include "auxiliary.h"
#include "encoding.h"
#include "lexmin_solver.h"
#include "options.h"
#include "read_gap.h"
#include "read_mace.h"
#include "statistics.h"
#include "trie.h"
#include "version.h"
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;
static void prn_header(Output &);
static void solve(Output &, const BinaryFunction &);
static void
solve_more(Output &                                            options,
           const std::vector<std::unique_ptr<BinaryFunction>> &tables);

int main(int argc, char **argv) {
    CLI::App          app("Lexicography smallest automorphic model.");
    Options           options;
    StatisticsManager statistics;
    Output            output(options, statistics);

    std::string file_name;
    app.add_option("file_name", file_name, "file name")->default_val("-");
    app.add_flag("-i", options.incremental, "Use incremental SAT solving.")
        ->default_val(true);
    app.add_flag("-v", options.verbose, "add verbosity")->default_val(0);
    app.add_flag("-u", options.unique, "output only unique models")
        ->default_val(0);
    app.add_flag("-m", options.mace_format, "use mace format for input/output")
        ->default_val(0);
    CLI11_PARSE(app, argc, argv);
    options.comment_prefix = options.mace_format ? "%" : "#";

    const bool use_std = file_name == "-";
    gzFile in = use_std ? gzdopen(0, "rb") : gzopen(file_name.c_str(), "rb");

    if (argc == 1)
        output.comment() << "Reading from standard input." << std::endl;

    output.comment(1) << "incrementality: "
                      << (options.incremental ? "true" : "false") << std::endl;
    output.comment(1) << "verbosity: " << options.verbose << std::endl;

    if (in == NULL) {
        printf("ERROR! Could not open file: %s\n",
               argc == 1 ? "<stdin>" : argv[1]);
        exit(EXIT_FAILURE);
    }
    prn_header(output);

    const auto start_time = read_cpu_time();
    if (options.mace_format) {
        ReadMace reader(output, in);
        reader.read();
        if (!use_std)
            gzclose(in);
        solve_more(output, reader.functions());
    } else {
        ReadGAP reader(in);
        reader.read();
        if (!reader.has_f()) {
            puts("function not read");
            exit(EXIT_FAILURE);
        }
        if (!use_std)
            gzclose(in);

        const auto &f = reader.f();
        solve(output, f);
    }
    statistics.totalTime->inc(read_cpu_time() - start_time);
    for (const auto s : statistics.all)
        s->print(output.comment()) << std::endl;
    return EXIT_SUCCESS;
}

static void
solve_more(Output &                                            output,
           const std::vector<std::unique_ptr<BinaryFunction>> &tables) {
    auto &                     options(output.d_options);
    std::unique_ptr<ModelTrie> mt(options.unique ? new ModelTrie() : nullptr);
    std::vector<std::unique_ptr<BinaryFunction>> unique_solutions;
    for (const auto &table : tables) {
        output.comment(1) << "solving " << table->get_name()
                          << " order:" << table->order() << " "
                          << table->get_additional_info() << std::endl;
        LexminSolver solver(output, *table);
        solver.solve();
        if (options.unique) {
            if (mt->add(*(solver.solution()))) {
                unique_solutions.push_back(std::move(solver.solution()));
                assert(solver.solution().get() == nullptr);
                output.comment(2)
                    << "added: " << unique_solutions.back()->order() << " "
                    << unique_solutions.back()->get_additional_info()
                    << std::endl;
            }
        } else {
            solver.print_solution(cout);
        }
    }
    if (options.unique) {
        for (const auto &table : unique_solutions)
            table->print_mace(cout);
    }
}

static void solve(Output &output, const BinaryFunction &table) {
    LexminSolver solver(output, table);
    solver.solve();
    solver.print_solution(cout);
}

void prn_header(Output &output) {
    output.comment() << "mlex, v00.0, " << Version::GIT_SHA1 << ", "
                     << Version::GIT_DATE << endl;
    output.comment() << "(C) 2022 Mikolas Janota, mikolas.janota@gmail.com"
                     << endl;
#ifdef USE_IPASIR
    output.comment() << "solver IPASIR (cadical)" << endl;
#endif /*USE_IPASIR*/
#ifdef USE_MINISAT
    output.comment() << "solver MINISAT" << endl;
#endif /*USE_MINISAT*/
#ifndef NDEBUG
    output.comment() << "DEBUG version (don't use for computationally heavy "
                        "instances)."
                     << endl;
#endif
}
