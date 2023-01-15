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
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;
static void prn_header(Output &);
static void solve(Output &, const BinaryFunction &);
static void solve_more(Output &options, ReadMace &reader);
static double start_time;

int main(int argc, char **argv) {
    CLI::App app("Lexicography smallest automorphic model.");
    Options options;
    StatisticsManager statistics;
    Output output(options, statistics);

    std::string file_name;
    app.add_option("file_name", file_name, "file name")->default_val("-");
    app.add_flag("-i", options.incremental, "Use incremental SAT solving.")
        ->default_val(true);
    app.add_flag("-1", options.opt1stRow,
                 "Try to optimize for the first row of the table.")
        ->default_val(false);
    app.add_flag("-b", options.budgeting, "budgeting.")->default_val(false);
    app.add_flag("-v", options.verbose, "add verbosity")->default_val(0);
    app.add_flag("-u", options.unique, "output only unique models")
        ->default_val(0);
    app.add_flag("-m", options.mace_format, "use mace format for input/output")
        ->default_val(0);
    app.add_flag("-r", options.invariants, "use row invariants")
        ->default_val(0);
    app.add_flag(
           "-l", options.last_solution,
           "Check last solution to see that this value is already possible.")
        ->default_val(false);
    app.add_option("--seq-counter-lits", options.seq_counter_lits,
                   "seq_counter_lits.")
        ->default_val(10);
    CLI11_PARSE(app, argc, argv);
    options.comment_prefix = options.mace_format ? "%" : "#";

    const bool use_std = file_name == "-";
    gzFile in = use_std ? gzdopen(0, "rb") : gzopen(file_name.c_str(), "rb");

    if (argc == 1)
        output.comment() << "Reading from standard input." << std::endl;

    output.comment(1) << "incrementality: "
                      << (options.incremental ? "true" : "false") << std::endl;
    output.comment(1) << "verbosity: " << options.verbose << std::endl;
    output.comment(1) << "seq_counter_lits: " << options.seq_counter_lits
                      << std::endl;

    if (in == nullptr) {
        printf("ERROR! Could not open file: %s\n",
               argc == 1 ? "<stdin>" : argv[1]);
        exit(EXIT_FAILURE);
    }
    prn_header(output);

    start_time = read_cpu_time();
    if (options.mace_format) {
        ReadMace reader(output, in);
        solve_more(output, reader);
        if (!use_std)
            gzclose(in);
    } else {
        ReadGAP reader(in);
        reader.read();
        statistics.readModels->inc();
        statistics.readingTime->inc(read_cpu_time() - start_time);
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
process_tables(Output &output,
               const std::vector<std::unique_ptr<BinaryFunction>> &tables,
               std::unique_ptr<ModelTrie> &mt,
               std::vector<std::unique_ptr<BinaryFunction>> &unique_solutions,
               size_t &counter) {
    auto &options(output.d_options);
    auto &statistics(output.d_statistics);

    for (const auto &table : tables) {
        statistics.readModels->inc();
        if (counter && (counter % 1000) == 0) {
            if (options.verbose == 0) {
                (output.ccomment()
                 << " " << counter << " ("
                 << SHOW_TIME0(read_cpu_time() - start_time) << "s)")
                    .flush();
            } else {
                output.comment() << "done: " << counter << std::endl;
            }
        }

        output.comment(1) << "solving " << table->get_name()
                          << " order:" << table->order() << " "
                          << table->get_additional_info() << std::endl;

        LexminSolver solver(output, *table);
        solver.solve();
        if (options.unique) {
            if (mt->add(*(solver.solution()))) {
                solver.solution()->set_name(table->get_name());
                solver.solution()->set_additional_info(
                    table->get_additional_info());
                unique_solutions.push_back(std::move(solver.solution()));
                assert(solver.solution().get() == nullptr);
                output.comment(2)
                    << "added: " << unique_solutions.back()->order() << " "
                    << unique_solutions.back()->get_additional_info()
                    << std::endl;
            }
        } else {
            statistics.producedModels->inc();
            solver.solution()->set_name(table->get_name());
            solver.solution()->set_additional_info(
                table->get_additional_info());
            solver.print_solution(cout);
        }
        counter++;
    }
}

static int read(Output &output, ReadMace &reader, int max_read) {
    const auto t0 = read_cpu_time();
    const auto rv = reader.read(max_read);
    output.d_statistics.readingTime->inc(read_cpu_time() - t0);
    return rv;
}

static void solve_more(Output &output, ReadMace &reader) {
    auto &options(output.d_options);
    auto &statistics(output.d_statistics);

    std::unique_ptr<ModelTrie> mt(options.unique ? new ModelTrie() : nullptr);
    std::vector<std::unique_ptr<BinaryFunction>> unique_solutions;

    size_t counter = 0;

    if (options.verbose == 0)
        output.comment() << "Done:";

    const auto max_read = 100;
    while (read(output, reader, max_read)) {
        process_tables(output, reader.functions(), mt, unique_solutions,
                       counter);
        reader.clear();
    };

    if (options.verbose == 0)
        output.ccomment() << std::endl;

    if (options.unique) {
        for (const auto &table : unique_solutions) {
            table->print_mace(cout);
            statistics.producedModels->inc();
        }
    }
}

static void solve(Output &output, const BinaryFunction &table) {
    LexminSolver solver(output, table);
    solver.solve();
    solver.print_solution(cout);
    output.d_statistics.producedModels->inc();
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
