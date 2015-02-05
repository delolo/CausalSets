/*----------------------------------------------------------------
 *
 *  Written:       13/01/2014
 *  Last updated:  13/01/2014
 *
 *
 *  This class contains the methods necessary to test Fay's
 *  conjecture: the integral of the extrinsic curvature over a
 *  hypersurface sigma is equal to the rho->infinity limit of
 *
 *  const * rho^{(2-d)/2} * (Min(rho,sigma) - Max(rho,sigma))
 *
 *  where:
 *  const depends on the spacetime dimension only;
 *  Min(rho, sigma) is the number of minimal causet elements to the
 *  future of sigma;
 *  Max(rho, sigma) is the number of maximal causet elements to the
 *  past of sigma;
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

class TExtrinsicCurvature {

private:
    int dimension;
    double totalvolume;
    TGeometry* GeometryAbove;
    TGeometry* GeometryBelow;
    double rhoexp;
    int expSizeAbove;
    int expSizeBelow;
    int currentSprinklingSize;

public:
    // constructor for fixed density analysis. specifies spacetime region below and above
    // the surface sigma.
    TExtrinsicCurvature(int dimension, TGeometry *GeometryAbove, TGeometry *GeometryBelow,
            double rhoexp) {
        this->dimension = dimension;
        this->GeometryAbove = GeometryAbove;
        this->GeometryBelow = GeometryBelow;
        this->totalvolume = GeometryAbove->V + GeometryBelow->V;
        this->rhoexp = rhoexp;
        this->expSizeAbove = (int) (rhoexp * (GeometryAbove->V));
        this->expSizeBelow = (int) (rhoexp * (GeometryBelow->V));
        this->currentSprinklingSize = 0;
    }

    // constructor for fixed number analysis. specifies spacetime region below and above
    // the surface sigma.
    TExtrinsicCurvature(int dimension, TGeometry *GeometryAbove, TGeometry *GeometryBelow,
            int nexp) {
        this->dimension = dimension;
        this->GeometryAbove = GeometryAbove;
        this->GeometryBelow = GeometryBelow;
        this->totalvolume = GeometryAbove->V + GeometryBelow->V;
        this->expSizeAbove = round(nexp / 2);
        this->expSizeBelow = round(nexp / 2);
        this->rhoexp = nexp / totalvolume;
        this->currentSprinklingSize = 0;
    }

    // empty destructor
    ~TExtrinsicCurvature() {
    }

    void setDensity(double density) {
        this->rhoexp = density;
        this->expSizeAbove = round(this->rhoexp * GeometryAbove->V / 2);
        this->expSizeBelow = round(this->rhoexp * GeometryBelow->V / 2);
    }

    // returns a vector with mean and standard deviation
    vector<double> getNewStats(int N) {
        vector<double> out;
        double mean = 0;
        double stdev = 0;
        double current;
        for (int i = 0; i < N; i++) {
            current = this->c(this->dimension)
                    * pow(rhoexp, 2.0 / dimension - 1.0)
                    * getNewMinMinusMax();
            mean += current;
            stdev += current * current;
        }
        mean = mean / N;
        stdev = sqrt((stdev - mean * mean * N) / (N - 1));
        out.push_back(mean);
        out.push_back(stdev);
        out.push_back(stdev / sqrt((double) N));
        return out;
    }

    // returns a vector with mean and standard deviation
    void printNewDataAndStats(int N, ofstream* data, ofstream* stat) {
        double mean = 0;
        double stdev = 0;
        double current;
        double meanelements = 0;
        std::clock_t start = std::clock();
        for (int i = 0; i < N; i++) {
            // print status to console
            if ((i + 1) % (int)(N / 100) == 0) {
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
            current = getNewMinMinusMax();
//                    * this->c(this->dimension)
//                    * pow(rhoexp, 2.0 / dimension - 1.0);
            meanelements += currentSprinklingSize;
            (*data) << this->rhoexp << "\t" << current << "\n";
            (*data).flush();
            mean += current;
            stdev += current * current;
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
    }
    ;

    // returns a vector with mean and standard deviation
    void printNewDataAndStatsDebug(int N, ofstream* data, ofstream* stat) {
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
            std::vector<int> minmax = getNewMinMinusMaxDebug();
            //current = minmax[0];
//                    * this->c(this->dimension)
//                    * pow(rhoexp, 2.0 / dimension - 1.0)
            meanelements += currentSprinklingSize;
            (*data) << this->rhoexp << "\t" << minmax[0] << "\t" << minmax[1] << "\n";
            (*data).flush();
            mean += minmax[0];
            stdev += minmax[0] * minmax[0];
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
    int getNewMinMinusMax() {
        TSprinkling* SprinklingAbove = new TSprinkling(expSizeAbove,
                GeometryAbove);
        TSprinkling* SprinklingBelow = new TSprinkling(expSizeBelow,
                GeometryBelow);
        this->currentSprinklingSize = SprinklingAbove->size
                + SprinklingBelow->size;
        int out = countMinimalElements(SprinklingAbove)
                - countMaximalElements(SprinklingBelow);
        delete SprinklingAbove;
        delete SprinklingBelow;
        return out;
    }

    std::vector<int> getNewMinMinusMaxDebug() {
        std::vector<int> out;
        TSprinkling* SprinklingAbove = new TSprinkling(expSizeAbove,
                GeometryAbove);
        TSprinkling* SprinklingBelow = new TSprinkling(expSizeBelow,
                GeometryBelow);
        this->currentSprinklingSize = SprinklingAbove->size
                + SprinklingBelow->size;
        out.push_back(countMinimalElements(SprinklingAbove));
        out.push_back(countMaximalElements(SprinklingBelow));
        delete SprinklingAbove;
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
            while (futureset < k && (j < Sprinkling->size)) {
                if (Sprinkling->prec(i, j)) futureset++;
                j++;
            }
            if (futureset == k) out++;
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
