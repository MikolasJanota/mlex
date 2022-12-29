/*
 * File:  statistics.h
 * Author:  mikolas
 * Created on:  Fri Dec 23 12:47:28 CET 2022
 * Copyright (C) 2022, Mikolas Janota
 */
#pragma once
#include <iostream>
#include <string>
#include <vector>

class StatisticsManager {
  public:
    StatisticsManager() {
        all.push_back(satCalls = new IntStatistic("SAT calls"));
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

      private:
        const std::string d_name;
    };

    class DoubleStatistic : public Statistic {
      public:
        DoubleStatistic(const std::string &name, double init_value = 0)
            : Statistic{name}, d_val{init_value} {};

        double inc(double ival) { return d_val += ival; }
        double get() const { return d_val; }

        virtual std::ostream &print(std::ostream &o) override {
            return o << name() << " : " << d_val;
        }

      private:
        double d_val;
    };

    class IntStatistic : public Statistic {
      public:
        IntStatistic(const std::string &name, int init_value = 0)
            : Statistic{name}, d_val{init_value} {};
        int inc() { return ++d_val; }
        int get() const { return d_val; }

        virtual std::ostream &print(std::ostream &o) override {
            return o << name() << " : " << d_val;
        }

      private:
        int d_val;
    };
    IntStatistic *           satCalls;
    DoubleStatistic *        satTime;
    DoubleStatistic *        readingTime;
    DoubleStatistic *        totalTime;
    std::vector<Statistic *> all;
};
