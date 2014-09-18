/*----------------------------------------------------------------
 *
 *  Written:       26/05/2014
 *  Last updated:  26/05/2014
 *
 *
 *  This class contains the methods necessary to test the smeared
 *  version of Fay's conjecture. See UnitExtrinsicCurvature for more
 *  info.
 *
 *  TODO: Current implementation has methods for 2d sprinklings only.
 *        Also it is restricted to t=const surfaces.
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
#include "UnitConnection.h"

using namespace std;

class TExtrinsicCurvatureSmeared {

private:
    int dimension;
    TGeometry* GeometryAbove;
    TGeometry* GeometryBelow;
    double rho;
    int expSizeAbove;
    int expSizeBelow;
    int currentSprinklingSize;

public:
    // constructor for fixed density analysis. specifies spacetime region below and above
    // the surface sigma.
    TExtrinsicCurvatureSmeared(int dimension, TGeometry *GeometryAbove, TGeometry *GeometryBelow, double density) {
        this->dimension = dimension;
        this->GeometryAbove = GeometryAbove;
        this->GeometryBelow = GeometryBelow;
        this->rho = density;
        this->expSizeAbove = (int) (rho * (GeometryAbove->V));
        this->expSizeBelow = (int) (rho * (GeometryBelow->V));
        this->currentSprinklingSize = 0;
    }

    // empty destructor
    ~TExtrinsicCurvatureSmeared() {
    }

    void setDensity(double density) {
        this->rho = density;
    }

    // returns a vector with mean and standard deviation
    vector<double> getNewStats(int N, double epsilon) {
        vector<double> out;
        double mean = 0;
        double stdev = 0;
        double current;
        for (int i = 0; i < N; i++) {
            current = pow(epsilon * rho, 2.0 / dimension - 1.0) * getNewMinMinusMax(epsilon);
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
    void printNewDataAndStats(int N, double epsilon, ofstream* data, ofstream* stat) {
        double mean = 0;
        double stdev = 0;
        double current;
        double meanelements = 0;
        std::clock_t start = std::clock();
        for (int i = 0; i < N; i++) {
            // print status to console
            if ((i + 1) % (N / 100) == 0) {
                cout << "Progress: "
                     << 100 * (double) (i + 1) / N << "\%"
                     << endl;
            }
            // print ETA to console
            if ((i + 1) == N / 20) {
                cout << "PROJECTED TIME = "
                     << 20 * (std::clock() - start) / (double) CLOCKS_PER_SEC
                     << "s" << endl;
            }
            // find stats
            current = pow(epsilon * rho, 2.0 / dimension - 1.0) * getNewMinMinusMax(epsilon);
            meanelements += currentSprinklingSize;
            (*data) << this->rho << "\t" << current << "\n";
            mean += current;
            stdev += current * current;
        }
        meanelements = meanelements / N;
        mean = mean / N;
        stdev = sqrt((stdev - mean * mean * N) / (N - 1));
        double totaltime = (std::clock() - start) / (double) CLOCKS_PER_SEC;
        (*stat) << std::setw(0) << N
                << std::setw(16) << std::fixed << std::setprecision(0) << rho
                << std::setw(16) << std::setprecision(3) << std::scientific << epsilon
                << std::setw(16) << std::fixed << meanelements
                << std::setw(10) << totaltime
                << std::setw(10) << std::fixed << std::setprecision(3) << mean
                << std::setw(10) << std::fixed << std::setprecision(3) << stdev
                << endl;
    }
    ;

public:
    double getNewMinMinusMax(double epsilon) {
        TSprinkling* SprinklingAbove = new TSprinkling(expSizeAbove, GeometryAbove);
        TSprinkling* SprinklingBelow = new TSprinkling(expSizeBelow, GeometryBelow);
        this->currentSprinklingSize = SprinklingAbove->size
                                    + SprinklingBelow->size;
        double out = countMinimalElements(SprinklingAbove, epsilon)
                   - countMaximalElements(SprinklingBelow, epsilon);
        delete SprinklingAbove;
        delete SprinklingBelow;
        return out;
    }

    double countMaximalElements(TSprinkling* Sprinkling, double epsilon) {
        if (Sprinkling->size == 0) return 0.0;
        TConnection* Connection = new TConnection(Sprinkling, false, false);
        /*Connection->writeCausalMatrix("c-below.dat");*/
        double out = 0.0;
        for(int i = 0; i < Connection->size; i++) {
            out += pow(1.0 - epsilon, Connection->cardOfPastSet(i));
        }
        delete Connection;
        return out;
    }

    double countMinimalElements(TSprinkling *Sprinkling, double epsilon) {
        if (Sprinkling->size == 0) return 0.0;
        TConnection* Connection = new TConnection(Sprinkling, false, false);
        /*Connection->writeCausalMatrix("c-above.dat");*/
        double out = 0.0;
        for(int i = 0; i < Connection->size; i++) {
            out += pow(1.0 - epsilon, Connection->cardOfFutureSet(i));
        }
        delete Connection;
        return out;
    }

};
