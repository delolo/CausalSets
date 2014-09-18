/*----------------------------------------------------------------
 *
 *  Written:       14/02/2014
 *  Last updated:  14/02/2014
 *
 *
 *  This is a statistics class. Use as follows:
 *
 *  Statistics myStats = new Statistics(...);
 *
 *  When fed with new data points x, the class keeps track of
 *  s1 = sum(x) and s2 = sum(x*x) (as described on e.g.
 *  Wikipedia) to calculate the running sample mean
 *
 *  mean = s1 / N
 *
 *  and the sample standard deviation
 *
 *  sdev = sqrt((N * s1 - s2 * s2)/(N * (N - 1))).
 *
 *----------------------------------------------------------------*/

#pragma once

#include <math.h>
#include "time.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "StringTools.h"

class Statistics {

private:
    /* boolean that specifies whether all data is stored in a local vector */
    bool keepdata;
    /* vector that optionally keeps all the data */
    std::vector<double> data;
    /* starting time */
    std::clock_t start;
    /* time of last call to clock */
    std::clock_t lastcall;
    /* number of calls made to clock */
    int calls;
    /* name of statistics object */
    std::string name;
    /* size[i] = number of data points for ith variable*/
    int s0;
    /* s1[i] = s1 for ith variable*/
    double s1;
    /* s1[i] = standard deviation for ith variable*/
    double s2;

public:
    /*  Constructor for one variable*/
    Statistics(const char* name) {
        this->calls = 0;
        this->s0 = 0;
        this->s1 = 0.0;
        this->s2 = 0.0;
        this->name = name;
        this->keepdata = false;
        this->data.resize(0);
        start = std::clock();
        lastcall = start;
        std::cout << "Starting statistics for \"" << name << "\"." << std::endl;
    }

    /*  Constructor for one variable*/
    Statistics(const char* name, bool _keepdata) {
        this->calls = 0;
        this->s0 = 0;
        this->s1 = 0.0;
        this->s2 = 0.0;
        this->name = name;
        this->keepdata = _keepdata;
        this->data.resize(0);
        start = std::clock();
        lastcall = start;
        std::cout << "Starting statistics for \"" << name << "\"." << std::endl;
    }



    /* Add new data-point to first variable */
    void push(double x) {
        if(keepdata) pushandkeep(x);
        else pushanddrop(x);
    }

private:
    void pushanddrop(double x) {
        this->s0++;
        this->s1 += x;
        this->s2 += x * x;
    }

    void pushandkeep(double x) {
        data.push_back(x);
        this->s0++;
        this->s1 += x;
        this->s2 += x * x;
    }

public:
    //  records time taken since instantiation and
    //  prints it out to the console
    std::string stringTimeTotal() {
        std::ostringstream stream;
        double t = (std::clock() - this->start) / (double) CLOCKS_PER_SEC;
        this->lastcall = std::clock();
        stream << "Process \"" << name << "\" - Total Time: " << t << "s."
                << std::endl;
        return stream.str();
    }

    //  records time taken since instantiation and
    //  prints it out to the console
    std::string getTime() {
        std::ostringstream stream;
        double t = (std::clock() - this->lastcall) / (double) CLOCKS_PER_SEC;
        this->lastcall = std::clock();
        this->calls++;
        stream  << std::setprecision(5)
                << "PROCESS \""
                << name
                << "\" - SPLIT "
                << calls
                << " TIME: "
                << t
                << "s."
                << std::endl;
        return stream.str();
    }

    /* get stats for first variable */
    std::vector<double> getStats() {
        std::vector<double> out(3);
        out[0] = s0;
        out[1] = s1 / s0;
        out[2] = sqrt((s0 * s2 - s1 * s1)/ (s0 * (s0 - 1)));
        return out;
    }

    double getMean() {
        return s1 / s0;
    }

    double getSdev() {
        return sqrt((s0 * s2 - s1 * s1)/ (s0 * (s0 - 1)));
    }

    /* toString() */
    std::string toString() {
        std::ostringstream stream;
        std::vector<double> results = this->getStats();
        stream << std::setw(8) << std::setprecision(3);
        stream << std::scientific;
        stream << (std::clock() - this->start) / (double) CLOCKS_PER_SEC << "\t";
        stream << std::fixed << std::setw(8) << std::setprecision(0) << std::right;
        stream << results[0] << "\t";
        stream << std::setw(8) << std::setprecision(3);
        stream << results[1] << "\t" << results[2] << "\n";
        return stream.str();
    }

    /* toStringFirstLine() */
    std::string toStringFirstLine() {
        std::ostringstream stream;
        stream << "TIME\tSAMPLE\tMEAN\tSDEV";
        return stream.str();
    }

    /* toStringFirstLine() */
    void dataToFile(const char* name) {
        if(!keepdata) throw "Data was not kept";
        std::ofstream os;
        os.open(name);
        for(int i = 0; i < data.size(); i++) {
            os << data[i] << "\n";
        }
    }

    ;
};


