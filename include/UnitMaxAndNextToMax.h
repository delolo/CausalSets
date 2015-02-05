/*
 * UnitMinAndNextToMin.h
 *
 *  Created on: Feb 3, 2015
 *      Author: michelbuck
 */

#ifndef UNITMINANDNEXTTOMIN_H_
#define UNITMINANDNEXTTOMIN_H_


/*----------------------------------------------------------------
 *
 *  Written:       13/01/2014
 *  Last updated:  13/01/2014
 *
 *
 *  This class contains the methods necessary to produce the number of
 *  minimal and next-to-minimal causet elements to the future of sigma
 *
 *----------------------------------------------------------------*/

#pragma once

#include <iterator>
#include <new>
#include <vector>
#include <math.h>
#include <ctime>
#include <iomanip>

#include "UnitGeometry.h"
#include "UnitSprinkling.h"
#include <stdexcept>

using namespace std;

class TMaxAndNextToMax {

private:
    int dimension;
    double totalvolume;
    TGeometry* GeometryBelow;
    double rhoexp;
    int expSizeBelow;
    int currentSprinklingSize;

public:
    // constructor for fixed density analysis. specifies spacetime region below and above
    // the surface sigma.
    TMaxAndNextToMax(int dimension, TGeometry *GeometryBelow, double rhoexp) {
        this->dimension = dimension;
        this->GeometryBelow = GeometryBelow;
        this->totalvolume = GeometryBelow->V;
        this->rhoexp = rhoexp;
        this->expSizeBelow = (int) (rhoexp * (GeometryBelow->V));
        this->currentSprinklingSize = 0;
    }

    // constructor for fixed number analysis. specifies spacetime region below and above
    // the surface sigma.
    TMaxAndNextToMax(int dimension, TGeometry *GeometryBelow,
            int nexp) {
        this->dimension = dimension;
        this->GeometryBelow = GeometryBelow;
        this->totalvolume = GeometryBelow->V;
        this->expSizeBelow = round(nexp);
        this->rhoexp = nexp / totalvolume;
        this->currentSprinklingSize = 0;
    }

    // empty destructor
    ~TMaxAndNextToMax() {
    }

    void setDensity(double density) {
        this->rhoexp = density;
        this->expSizeBelow = round(this->rhoexp * GeometryBelow->V);
    }

    void printNewMaxAndNextToMaxDataAndStats(int N, ofstream* data, ofstream* stat) {
        double mean = 0;
        double stdev = 0;
        double meanelements = 0;
        std::clock_t start = std::clock();
        for (int i = 0; i < N; i++) {
            // print status to console
            if ((i + 1) % (N / 100) == 0) {
                cout << "Progress: "
                        << 100 * (double) (i + 1) / N
                        << "\%"
                        << endl;
            }
            // print ETA to console
            if ((i + 1) == N / 20) {
                cout << "PROJECTED TIME = "
                        << 20 * (std::clock() - start) / (double) CLOCKS_PER_SEC
                        << "s" << endl;
            }
            // find stats
            std::vector<int> maxandnexttomax = getNewMaxAndNextToMax();
            //current = minmax[0];
//                    * this->c(this->dimension)
//                    * pow(rhoexp, 2.0 / dimension - 1.0)
            meanelements += currentSprinklingSize;
            (*data) << this->rhoexp << "\t" << maxandnexttomax[0] << "\t" << maxandnexttomax[1] << "\n";
            (*data).flush();
            mean += maxandnexttomax[0];
            stdev += maxandnexttomax[0] * maxandnexttomax[0];
        }
        meanelements = meanelements / N;
        mean = mean / N;
        stdev = sqrt((stdev - mean * mean * N) / (N - 1));
        double totaltime = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        (*stat) << std::setw(0)
                << N
                << std::setw(16) << std::fixed << std::setprecision(4)
                << this->rhoexp
                << std::setw(16) << std::fixed << std::setprecision(4)
                << this->rhoexp * this->totalvolume
                << std::setw(16) << std::fixed << std::setprecision(3)
                << meanelements
                << std::setw(10) << std::fixed << std::setprecision(3)
                << totaltime
                << std::setw(10) << std::fixed << std::setprecision(3)
                << mean
                << std::setw(10) << std::fixed << std::setprecision(3)
                << stdev
                << endl;
    };

private:
    std::vector<int> getNewMaxAndNextToMax() {
        TSprinkling* SprinklingBelow = new TSprinkling(expSizeBelow, GeometryBelow);
        std::vector<int> out = countMaximalAndSubMaximalElements(SprinklingBelow);
        delete SprinklingBelow;
        return out;
    }

    int countMaximalElements(TSprinkling* Sprinkling) {
        if (Sprinkling->size == 0) return 0;
        if (Sprinkling->size == 1) return 1;
        int out = 1; // the last element is always maximal
        bool ismax;
        int j;
        for (int i = 0; i < Sprinkling->size - 1; i++) {
            ismax = true;
            j = i + 1;
            while (ismax && (j < Sprinkling->size)) {
                if (Sprinkling->prec(i, j)) ismax = false;
                j++;
            }
            if (ismax) out++;
        }
        return out;
    }

    int countSubMaximalElements(TSprinkling* Sprinkling, int k) {
        if (Sprinkling->size == 0) return 0;
        if (k == 0) return countMaximalElements(Sprinkling);
        int out = 0;
        int j = 0;
        int futureset; // counts the elements to the future of i
        /*  Elements n-k,...,n cannot be k-to-maximal
         *  so only iterate up to n - k - 1           */
        for (int i = 0; i < Sprinkling->size - k; i++) {
            futureset = 0;
            j = i + 1;
            while (futureset <= k && (j < Sprinkling->size)) {
                if (Sprinkling->prec(i, j)) futureset++;
                j++;
            }
            if (futureset == k) out++;
        }
        return out;
    }

    std::vector<int> countMaximalAndSubMaximalElements(TSprinkling* Sprinkling) {
            std::vector<int> out;
            out.push_back(0);
            out.push_back(0);
            int j = 0;
            int futureset; // counts the elements to the future of i
            for (int i = 0; i < Sprinkling->size; i++) {
                futureset = 0;
                j = i + 1;
                while (futureset <= 1 && (j < Sprinkling->size)) {
                    if (Sprinkling->prec(i, j)) futureset++;
                    j++;
                }
                if (futureset == 0) out[0]++;
                else if (futureset == 1) out[1]++;
            }
            return out;
        }

    int countMinimalElements(TSprinkling *Sprinkling) {
        if (Sprinkling->size == 0) return 0;
        if (Sprinkling->size == 1) return 1;
        int out = 1; // the first element is always minimal
        bool ismin;
        int j;
        for (int i = 1; i < Sprinkling->size; i++) {
            ismin = true;
            j = 0;
            while (ismin && (j < i)) {
                if (Sprinkling->prec(j, i)) ismin = false;
                j++;
            }
            if (ismin) out++;
        }
        return out;
    }

    int countSubMinimalElements(TSprinkling* Sprinkling, int k) {
        if (Sprinkling->size == 0) return 0;
        if (k == 0) return countMaximalElements(Sprinkling);
        int out = 0;
        int j;
        int pastset; // counts the elements to the future of i
        /*  Elements 0,...,k-1 cannot be k-to-minimal
         *  so only iterate from i = k */
        for (int i = k; i < Sprinkling->size; i++) {
            pastset = 0;
            j = 0;
            while (pastset < k && (j < i)) {
                if (Sprinkling->prec(j, i)) pastset++;
                j++;
            }
            if (pastset == k) out++;
        }
        return out;
    }

    double c(int dim) {
        if (dim == 2) return 0.75000;
        if (dim == 3) return 0.91386;
        if (dim == 4) return 0.96225;
        if (dim == 5) return 0.96099;
        else throw std::invalid_argument("I only know c(d) for d = 2 to d = 5.");
    }
};



#endif /* UNITMINANDNEXTTOMIN_H_ */
