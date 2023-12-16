/*
 * File:  statistics.h
 * Author:  mikolas
 * Created on:  Fri Dec 23 12:47:28 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class StatisticsManager {
  public:
    StatisticsManager() {
        all.push_back(producedModels = new IntStatistic("produced models"));
        all.push_back(readModels = new IntStatistic("read models"));
        all.push_back(unique1stRow = new IntStatistic("unique 1st row"));
        all.push_back(unique1stRowFirst =
                          new IntStatistic("unique 1st row is the first row"));
        all.push_back(fixedElements = new IntStatistic("fixed elements"));
        all.push_back(disallowed = new IntStatistic("disallowed mappings"));
        all.push_back(uniqueInv = new IntStatistic("unique by invariants"));
        all.push_back(maxCol = new IntStatistic("maximum color"));
        all.push_back(uniqueColInv =
                          new IntStatistic("unique by col invariants"));
        all.push_back(uniqueDiag1Elem =
                          new IntStatistic("unique diag 1st elem"));
        all.push_back(uniqueDiagElem = new IntStatistic("unique diag elem"));
        all.push_back(uniqueRowElem = new IntStatistic("unique row elem"));
        all.push_back(inferredCells = new IntStatistic("inferred cells"));
        all.push_back(satCalls = new IntStatistic("SAT calls"));
        all.push_back(encodingTime = new DoubleStatistic("encoding time"));
        all.push_back(satTime = new DoubleStatistic("SAT time"));
        all.push_back(readingTime = new DoubleStatistic("Reading time"));
        all.push_back(totalTime = new DoubleStatistic("Total time"));
    }

    virtual ~StatisticsManager();
    class Statistic {
      public:
        Statistic(const std::string &name);
        virtual ~Statistic();

        const std::string &name() const { return d_name; }

        virtual std::ostream &print(std::ostream &) = 0;
        virtual bool should_print() const = 0;

      private:
        const std::string d_name;
    };

    class DoubleStatistic : public Statistic {
      public:
        DoubleStatistic(const std::string &name, double init_value = 0)
            : Statistic{name}, d_val{init_value} {};

        double inc(double ival) {
            d_print = true;
            return d_val += ival;
        }

        double get() const { return d_val; }

        virtual bool should_print() const override { return d_print; }

        virtual std::ostream &print(std::ostream &o) override {
            return o << name() << " : " << std::fixed << std::setprecision(3)
                     << d_val;
        }

      private:
        double d_val;
        bool d_print = false;
    };

    class IntStatistic : public Statistic {
      public:
        IntStatistic(const std::string &name, int init_value = 0)
            : Statistic{name}, d_val{init_value} {};
        int inc() { return ++d_val; }
        int get() const { return d_val; }
        void set(int v) { d_val = v; }

        virtual bool should_print() const override { return d_val != 0; }

        virtual std::ostream &print(std::ostream &o) override {
            return o << name() << " : " << d_val;
        }

      private:
        int d_val;
    };
    IntStatistic *satCalls;
    IntStatistic *unique1stRow;
    IntStatistic *unique1stRowFirst;
    IntStatistic *fixedElements;
    IntStatistic *disallowed;
    IntStatistic *uniqueDiagElem;
    IntStatistic *uniqueRowElem;
    IntStatistic *inferredCells;
    IntStatistic *uniqueDiag1Elem;
    IntStatistic *uniqueInv;
    IntStatistic *uniqueColInv;
    IntStatistic *maxCol;
    IntStatistic *producedModels;
    IntStatistic *readModels;
    DoubleStatistic *encodingTime;
    DoubleStatistic *satTime;
    DoubleStatistic *readingTime;
    DoubleStatistic *totalTime;
    std::vector<Statistic *> all;
};
