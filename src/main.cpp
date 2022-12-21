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

static void
solve_more(const Options &                                     options,
           const std::vector<std::unique_ptr<BinaryFunction>> &tables) {
    for (const auto &table : tables) {
        LexminSolver solver(options, *table);
        solver.solve();
        solver.print_solution(cout);
    }
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
