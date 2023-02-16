/*
 * File:  main.cpp
 * Author:  mikolas
 * Created on:  Tue Dec 13 12:04:23 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#include "CLI11.hpp"
#include "auxiliary.h"
#include "binary_function.h" // for BinaryFunction
#include "comp_function.h"   // for CompFunction, CompFunction_hash, CompFu...
#include "graph.h"
#include "lexmin_solver.h"
#include "options.h"
#include "read_diags.h"
#include "read_gap.h"
#include "read_mace.h"
#include "statistics.h"
#include "trie.h"
#include "version.h"
#include <cassert> // for assert
#include <cstddef>
#include <cstdint> // for uint64_t
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>       // for map
#include <memory>    // for unique_ptr
#include <stdexcept> // for invalid_argument, out_of_range
#include <string>
#include <type_traits> // for remove_reference<>::type
#include <utility>     // for move
#include <vector>      // for vector
#include <zlib.h>      // for gzclose, gzdopen, gzopen, gzFile

using namespace std;
static void prn_header(Output &);
static void solve_more(Output &options, ReadMace &reader);
static void prn_mace(Output &output, ReadMace &reader);
static void solve_more_gaps(Output &options, ReadGAP &reader, ReadDiags *);
static double start_time;
/* #define LM_SET */

#ifndef LM_SET
#include <unordered_set>
using HT =
    std::unordered_set<CompFunction, CompFunction_hash, CompFunction_equal>;
#else
#include <set>
using HT = std::set<CompFunction, CompFunction_less>;
#endif

int main(int argc, char **argv) {
    CLI::App app("Lexicography smallest automorphic model.");
    Options options;
    StatisticsManager statistics;
    Output output(options, statistics);

    app.add_option("file_name", options.file_name, "file name")
        ->default_val("-");
    app.add_flag("-1, !--no-1", options.opt1stRow,
                 "Try to optimize for the first row of the table")
        ->default_val(true);
    app.add_flag("-b,!--no-b", options.budgeting, "Budgeting")
        ->default_val(true);
    app.add_flag("-v", options.verbose, "Add verbosity")->default_val(0);
    app.add_option("--diag-file", options.diag_file, "diag_file")
        ->default_val("");
    app.add_flag("-G", options.graph, "Graph")->default_val(false);
    app.add_flag("-P", options.print, "Print into files")->default_val(false);
    app.add_flag("-m", options.mace_format, "Use mace format for input/output")
        ->default_val(false);
    app.add_flag("-u", options.unique, "Output only unique models")
        ->default_val(false);

    // adapted from
    // https://github.com/CLIUtils/CLI11/blob/main/examples/enum.cpp
    std::map<std::string, SearchType> map{{"lus", SearchType::lin_us},
                                          {"lsu", SearchType::lin_su},
                                          {"bin", SearchType::bin},
                                          {"bin2", SearchType::bin2}};
    // CheckedTransformer translates and checks whether the results are either
    // in one of the strings or in one of the translations already
    app.add_option("-t,--search-type", options.search_type,
                   "set the search type")
        ->default_val(SearchType::bin2)
        ->transform(CLI::CheckedTransformer(map, CLI::ignore_case));
    ///// TODO better help message

    app.add_flag("-d,!--no-d", options.diagonal, "Traverse diagonal first")
        ->default_val(false);
    app.add_flag("-r,!--no-r", options.invariants, "Use row invariants")
        ->default_val(true);
    app.add_flag("-H,!--no-H", options.use_hash_table,
                 "Use hashtable to store unique models instead of trie")
        ->default_val(true);
    app.add_flag(
           "-l,!--no-l", options.last_solution,
           "Check last solution to see that this value is already possible")
        ->default_val(true);
    app.add_option(
           "--seq-counter-lits", options.seq_counter_lits,
           "Number of literals when to switch to seq counter enc for at most 1")
        ->default_val(10);
    app.add_flag("--simp_sat_row,!--no-simp_sat_row", options.simp_sat_row,
                 "Simplify SAT solver each row (minisat)")
        ->default_val(false);

    CLI11_PARSE(app, argc, argv);
    options.comment_prefix = options.mace_format ? "%" : "#";

    const bool use_std = options.file_name == "-";
    gzFile in =
        use_std ? gzdopen(0, "rb") : gzopen(options.file_name.c_str(), "rb");

    if (use_std)
        output.comment() << "Reading from standard input." << endl;
    else
        output.comment() << "Reading from " << options.file_name << endl;

    output.comment(1) << "verbosity: " << options.verbose << endl;
    output.comment(1) << "seq_counter_lits: " << options.seq_counter_lits
                      << std::endl;

    if (in == nullptr) {
        cerr << "ERROR! Could not open file: " << options.file_name << endl;
        exit(EXIT_FAILURE);
    }

    if (options.search_type == SearchType::lin_su && !options.last_solution) {
        cerr << "ERROR! multishot search requires last solution to be set"
             << endl;
        ;
        exit(EXIT_FAILURE);
    }
    if (!options.diag_file.empty() && !options.diagonal) {
        cerr << "ERROR! diagonal file only makes sense in the diagonal mode"
             << endl;
        ;
        exit(EXIT_FAILURE);
    }

    prn_header(output);

    start_time = read_cpu_time();

    if (options.mace_format) {
        ReadMace reader(output, in);
        if (options.graph || options.print) {
            prn_mace(output, reader);
        } else {
            solve_more(output, reader);
        }
        if (!use_std)
            gzclose(in);
    } else {
        gzFile din;
        ReadGAP reader(output, in);
        unique_ptr<ReadDiags> drd;
        if (!options.diag_file.empty()) {
            din = gzopen(options.diag_file.c_str(), "rb");
            drd = make_unique<ReadDiags>(output, din);
        }
        solve_more_gaps(output, reader, drd.get());
        gzclose(in);
        if (drd)
            gzclose(din);
    }
    statistics.totalTime->inc(read_cpu_time() - start_time);
    for (const auto s : statistics.all)
        s->print(output.comment()) << std::endl;
    return EXIT_SUCCESS;
}

static void process_table_plain(Output &output, const BinaryFunction &table,
                                LexminSolver &solver) {
    auto &statistics(output.d_statistics);
    const auto &options(output.d_options);
    statistics.producedModels->inc();
    std::unique_ptr<BinaryFunction> solution(solver.make_solution());
    solution->set_name(table.get_name());
    solution->set_additional_info(table.get_additional_info());
    if (options.mace_format)
        solution->print_mace(cout) << endl;
    else
        solution->print_gap(cout) << endl;
}

static void process_table_trie(
    Output &output, const BinaryFunction &table, LexminSolver &solver,
    std::unique_ptr<ModelTrie> &mt,
    std::vector<std::unique_ptr<BinaryFunction>> &unique_solutions) {
    std::unique_ptr<BinaryFunction> solution(solver.make_solution());
    if (mt->add(*solution)) {
        solution->set_name(table.get_name());
        solution->set_additional_info(table.get_additional_info());
        unique_solutions.push_back(std::move(solution));
        assert(solution.get() == nullptr);
        output.comment(2) << "added: " << unique_solutions.back()->order()
                          << " "
                          << unique_solutions.back()->get_additional_info()
                          << std::endl;
    }
}

static void process_table_ht(Output &output, const BinaryFunction &table,
                             LexminSolver &solver, std::unique_ptr<HT> &ht) {
    auto &statistics(output.d_statistics);
    const auto &options(output.d_options);
    CompFunction sol = solver.make_solution_comp();
    sol.set_name(table.get_name());
    /* sol.set_additional_info(table.get_additional_info()); */
    auto [it, successful] = ht->insert(sol);
    if (!successful) {
        sol.free();
        return;
    }
    if (options.mace_format)
        sol.print_mace(cout, table.get_additional_info()) << endl;
    else
        sol.print_gap(cout) << endl;
    statistics.producedModels->inc();
}

static void
process_tables(Output &output,
               const std::vector<std::unique_ptr<BinaryFunction>> &tables,
               std::unique_ptr<HT> &ht, std::unique_ptr<ModelTrie> &mt,
               std::vector<std::unique_ptr<BinaryFunction>> &unique_solutions,
               size_t &counter) {
    auto &options(output.d_options);
    auto &statistics(output.d_statistics);
    const auto use_ht = options.unique && options.use_hash_table;
    const auto use_trie = options.unique && !options.use_hash_table;

    for (const auto &table : tables) {
        statistics.readModels->inc();
        /* if (counter && (counter % 1000) == 0) { */
        /*     if (options.verbose == 0) { */
        /*         (output.ccomment() */
        /*          << " " << counter << " (" */
        /*          << SHOW_TIME0(read_cpu_time() - start_time) << "s)") */
        /*             .flush(); */
        /*     } else { */
        /*         output.comment() << "done: " << counter << std::endl; */
        /*     } */
        /* } */

        output.comment(1) << "solving " << table->get_name()
                          << " order:" << table->order() << " "
                          << table->get_additional_info() << std::endl;

        LexminSolver solver(output, *table);
        solver.solve();
        if (use_ht) {
            process_table_ht(output, *table, solver, ht);
        } else if (use_trie) {
            process_table_trie(output, *table, solver, mt, unique_solutions);
        } else {
            process_table_plain(output, *table, solver);
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

static int read(Output &output, ReadGAP &reader, int max_read) {
    const auto t0 = read_cpu_time();
    const auto rv = reader.read(max_read);
    output.d_statistics.readingTime->inc(read_cpu_time() - t0);
    return rv;
}

static void prn_mace(Output &output, ReadMace &reader) {
    auto &options(output.d_options);
    auto &statistics(output.d_statistics);
    for (size_t cnt = 0; read(output, reader, 1) > 0; ++cnt) {
        statistics.readModels->inc();
        if (reader.functions().empty()) {
            cerr << "Function not read" << endl;
            exit(EXIT_FAILURE);
        }
        const auto &f = *(reader.functions().begin()->get());

        const string output_file_name = options.file_name + "_" +
                                        std::to_string(cnt) +
                                        (options.graph ? ".dre" : ".gap");
        ofstream output_file(output_file_name);
        if (options.graph) {
            Graph g(f);
            g.make();
            g.print_nauty(output_file) << endl;
            /* g.print_dot(output_file) << endl; */
        } else {
            f.print_mace(output_file) << endl;
            /* if (dgs) { */
            /*     const string diag_output_file_name = */
            /*         options.file_name + "_" + std::to_string(cnt) + ".diag";
             */
            /*     ofstream diag_output_file(diag_output_file_name); */
            /*     print_vec(diag_output_file << "[[", min_diag, 1) */
            /*         << ", ()]]\n"; */
            /* } */
        }
        reader.clear();
    }
}

static void solve_more_gaps(Output &output, ReadGAP &reader, ReadDiags *dgs) {
    auto &options(output.d_options);
    auto &statistics(output.d_statistics);
    const auto use_ht = options.unique && options.use_hash_table;
    const auto use_trie = options.unique && !options.use_hash_table;
    std::unique_ptr<ModelTrie> mt(use_trie ? new ModelTrie() : nullptr);
    std::unique_ptr<HT> ht(use_ht ? new HT() : nullptr);
    std::vector<std::unique_ptr<BinaryFunction>> unique_solutions;

    for (size_t cnt = 0; read(output, reader, 1) > 0; ++cnt) {
        statistics.readModels->inc();
        if (reader.functions().empty()) {
            cerr << "Function not read" << endl;
            exit(EXIT_FAILURE);
        }
        const auto &f = *(reader.functions().begin()->get());
        std::vector<size_t> min_diag;
        if (dgs) {
            dgs->read(1);
            if (dgs->diags().empty()) {
                cerr << "diags not read" << endl;
                exit(EXIT_FAILURE);
            }
            min_diag = std::move(dgs->diags()[0]);
            if (min_diag.size() != f.order()) {
                cerr << "diagonal order mismatch (" << min_diag.size() << ")"
                     << endl;
                exit(EXIT_FAILURE);
            }
            dgs->clear();
        }

        if (options.graph || options.print) {
            const string output_file_name = options.file_name + "_" +
                                            std::to_string(cnt) +
                                            (options.graph ? ".dre" : ".gap");
            ofstream output_file(output_file_name);
            if (options.graph) {
                Graph g(f);
                g.make();
                g.print_nauty(output_file) << endl << endl;
                /* g.print_dot(output_file) << endl; */
            } else {
                f.print_gap(output_file << '[') << ']' << endl;
                if (dgs) {
                    const string diag_output_file_name =
                        options.file_name + "_" + std::to_string(cnt) + ".diag";
                    ofstream diag_output_file(diag_output_file_name);
                    print_vec(diag_output_file << "[[", min_diag, 1)
                        << ", ()]]\n";
                }
            }
            reader.clear();
            continue;
        }

        if (options.verbose > 3)
            f.print(cout, options.comment_prefix);
        LexminSolver solver(output, f);

        if (dgs) {
            solver.set_diag(min_diag);
        }
        output.comment(1) << "solving order:" << f.order() << " " << std::endl;

        solver.solve();

        if (use_ht) {
            process_table_ht(output, f, solver, ht);
        } else if (use_trie) {
            process_table_trie(output, f, solver, mt, unique_solutions);
        } else {
            process_table_plain(output, f, solver);
        }
        reader.clear();
        cnt++;
    }
}

static void solve_more(Output &output, ReadMace &reader) {
    auto &options(output.d_options);
    auto &statistics(output.d_statistics);
    const auto use_ht = options.unique && options.use_hash_table;
    const auto use_trie = options.unique && !options.use_hash_table;
    std::unique_ptr<ModelTrie> mt(use_trie ? new ModelTrie() : nullptr);
    std::unique_ptr<HT> ht(use_ht ? new HT() : nullptr);
    std::vector<std::unique_ptr<BinaryFunction>> unique_solutions;

    size_t counter = 0;

    /* if (options.verbose == 0) */
    /*     output.comment() << "Done:"; */

    const auto max_read = 100;
    while (read(output, reader, max_read)) {
        process_tables(output, reader.functions(), ht, mt, unique_solutions,
                       counter);
        reader.clear();
    }

    /* if (options.verbose == 0) */
    /*     output.ccomment() << std::endl; */

    if (use_ht) {
        /* for (const CompFunction &table : *ht) { */
        /*     table.print_mace(cout); */
        /*     statistics.producedModels->inc(); */
        /* } */
#ifndef NDEBUG
        std::vector<uint64_t *> cleanup;
        for (const CompFunction &table : *ht)
            cleanup.push_back(table.get_data());
        for (auto d : cleanup)
            delete[] d;
#endif
    }

    if (use_trie) {
        for (const auto &table : unique_solutions) {
            table->print_mace(cout);
            statistics.producedModels->inc();
        }
    }
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
