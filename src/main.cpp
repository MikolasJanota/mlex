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
#include "trie.h"
#include "version.h"
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;
static void prn_header();
static void solve(const Options &, const BinaryFunction &);
static void
solve_more(const Options &                                     options,
           const std::vector<std::unique_ptr<BinaryFunction>> &tables);

int main(int argc, char **argv) {
    prn_header();
    CLI::App    app("Lexicography smallest automorphic model.");
    Options     options;
    std::string file_name;
    app.add_option("file_name", file_name, "file name")->default_val("-");
    app.add_flag("-i", options.incremental, "Use incremental SAT solving.")
        ->default_val(0);
    app.add_flag("-v", options.verbose, "add verbosity")->default_val(1);
    app.add_flag("-u", options.unique, "output only unique models")
        ->default_val(0);
    app.add_flag("-m", options.mace_format, "use mace format for input/output")
        ->default_val(0);
    CLI11_PARSE(app, argc, argv);

    const bool use_std = file_name == "-";
    gzFile in = use_std ? gzdopen(0, "rb") : gzopen(file_name.c_str(), "rb");
    if (argc == 1)
        puts("c Reading from standard input.");

    if (in == NULL) {
        printf("ERROR! Could not open file: %s\n",
               argc == 1 ? "<stdin>" : argv[1]);
        exit(EXIT_FAILURE);
    }

    if (options.mace_format) {
        ReadMace reader(in);
        reader.read();
        if (!use_std)
            gzclose(in);
        solve_more(options, reader.functions());
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
        solve(options, f);
    }
    return EXIT_SUCCESS;
}

std::ostream &comment(const Options &options) {
    return std::cout << (options.mace_format ? '%' : '#') << " ";
}

static void
solve_more(const Options &                                     options,
           const std::vector<std::unique_ptr<BinaryFunction>> &tables) {
    std::unique_ptr<ModelTrie> mt(options.unique ? new ModelTrie() : nullptr);
    std::vector<std::unique_ptr<BinaryFunction>> unique_solutions;
    const auto                                   start_time = read_cpu_time();
    for (const auto &table : tables) {
        LexminSolver solver(options, *table);
        solver.solve();
        if (options.unique) {
            if (mt->add(*(solver.solution()))) {
                unique_solutions.push_back(std::move(solver.solution()));
                assert(solver.solution().get() == nullptr);
                if (options.verbose)
                    comment(options)
                        << "added: " << unique_solutions.back()->order() << " "
                        << unique_solutions.back()->get_additional_info()
                        << std::endl;
            }
        } else {
            solver.print_solution(cout);
        }
    }
    if (options.unique)
        for (const auto &table : unique_solutions)
            table->print_mace(cout);
    std::cout << " Total time : " << SHOW_TIME(read_cpu_time() - start_time)
              << std::endl;
}

static void solve(const Options &options, const BinaryFunction &table) {
    LexminSolver solver(options, table);
    solver.solve();
    solver.print_solution(cout);
}

void prn_header() {
    cout << "# mlex, v00.0, " << Version::GIT_SHA1 << ", " << Version::GIT_DATE
         << endl;
    cout << "# (C) 2022 Mikolas Janota, mikolas.janota@gmail.com" << endl;
#ifdef USE_IPASIR
    cout << "# solver IPASIR (cadical)" << endl;
#endif /*USE_IPASIR*/
#ifdef USE_MINISAT
    cout << "# solver MINISAT" << endl;
#endif /*USE_MINISAT*/
#ifndef NDEBUG
    cout << "# DEBUG version (don't use for computationally heavy "
            "instances)."
         << endl;
#endif
}
