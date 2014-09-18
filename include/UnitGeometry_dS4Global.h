/*************************************************************************
 *
 *  The UnitGeometry class file for a finite region of 3-dimensional dS in
 *  closed global coordinates with conformal time.
 *
 *  The volume form in these coordinates is
 *
 *      alpha^4 sec(t)^4 x 3-sphere-element
 *
 *  where alpha is the dS-radius. It's hard to generate a random variable
 *  with PDF sec(t)^4 so we use instead the global time coordinate tau
 *  defined by tan(t/2) = tanh(t/2 alpha) in which the t-factor of the volume
 *  element is cosh(t)^3. This is well approximated by exp(3 * tau) / 8 for
 *  large tau (in which most of the points lie as eta0 approaches PI/2).
 *
 *************************************************************************/

#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Tools.h"

class Tp_dS4Global: public TPoint {
public:
    int dimension;
    std::vector<double> x;

    Tp_dS4Global(double _t, vector<double> _x) {
        dimension = 4;
        t = _t;
        x = _x;
    }

    Tp_dS4Global() {
        dimension = 4;
        t = 0.0;
        x.resize(dimension);
        fill(x.begin(), x.end(), 0.0);
    }

    void write(ofstream *os) {
        (*os) << t << "\t";
        for (int i = 0; i < dimension; i++) {
            (*os) << x[i] << "\t";
        }
        (*os) << endl;
    }
};

//  The geometry class for a slab (0,t0) of dS4 in closed global coordinates
class TGeometry_dS4Global: public TGeometry {
public:
    int dimension;
    double t0, alpha;  // parameters for the geometry
    double tau0; // useful constants for getRandomPoint()

    // Constructor for target average degree k, cardinality n, and time interval [0,t0]
    TGeometry_dS4Global(double _t0) {
        dimension = 4;
        alpha = 1.0;
        t0 = _t0;
        tau0 = 2.0 * alpha * atanh(tan(t0 / 2.0));
        V = 2.0 / 3.0 * PI * PI * tan(t0) * (2.0 + pow(cos(t0), -2.0));
        cout << "V = " << V << endl;
    }

    // Get a random point according to the volume form.
    TPoint* getRandomPoint() {
        double tau = log((exp(3.0 * tau0) - 1.0) * rnd() + 1.0) / 3.0;
        double t = 2.0 * atan(tanh(tau / 2.0 / alpha));
        vector<double> x = rndSpherical(dimension, alpha);
        return new Tp_dS4Global(t, x);
    }

    // causal relations
    bool prec(TPoint *a, TPoint *b) {
        return (b->t - a->t > arcDist(a, b));
    }

    // spatial separation (length of the great circle segment)
    // between two points in ds3
    double arcDist(TPoint* _a, TPoint* _b) {
        Tp_dS4Global *a = (Tp_dS4Global*) _a;
        Tp_dS4Global *b = (Tp_dS4Global*) _b;

        double sigma = 0.0;
        for(int i = 0; i < 4; i++){
          sigma += a->x[i] * b->x[i];
        }
        return acos(sigma);
    }

};

