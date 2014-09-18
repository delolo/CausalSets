/*************************************************************************
 *
 *  The UnitGeometry class file for a EinsteinStaticUniverse in D spacetime dimensions
 *
 *************************************************************************/

#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Tools.h"

class Tp_EinsteinStaticUniverseD: public TPoint {
public:
    int dimension;
    double t;
    std::vector<double> x;

    Tp_EinsteinStaticUniverseD(int _dimension, double _t, vector<double> _x) {
        if (_x.size() != _dimension) throw "size of vector and dimension do not agree";
        dimension = _dimension;
        t = _t;
        x = _x;
    }
    Tp_EinsteinStaticUniverseD(double _dimension) {
        t = 0.0;
        dimension = _dimension;
        x.resize(dimension);
        x[0] = 1.0;
        for (unsigned i = 1; i < dimension; i++) {
            x[i] = 0.0;
        }
    }
    void write(ofstream *os) {
        (*os) << t << "\t";
        for (int i = 0; i < dimension; i++) {
            (*os) << x[i] << "\t";
        }
        (*os) << endl;
    }
};

//  A geometry class for sprinkling into a region of
//  the universal cover of ads2 in conformally flat coordinates
//  for times between 0 and 1
class TGeometry_EinsteinStaticUniverseD: public TGeometry {
public:
    int dimension;
    double tmin, tmax, radius;

    // constructor with radius and time interval [0,1]
    TGeometry_EinsteinStaticUniverseD(int _dimension, double _radius,
            double _tmin, double _tmax) {
        dimension = _dimension;
        radius = _radius;
        tmin = _tmin;
        tmax = _tmax;
        V = (tmax - tmin) * volSphere(dimension - 1, radius);
    }
// constructor with radius but time interval [0,1]
    TGeometry_EinsteinStaticUniverseD(int _dimension, double _radius) {
        dimension = _dimension;
        radius = _radius;
        tmin = 0.0;
        tmax = 1.0;
        V = (tmax - tmin) * volSphere(dimension - 1, radius);
    }
    // constructor with radius = 1 and time interval [0,1]
    TGeometry_EinsteinStaticUniverseD(int _dimension) {
        dimension = _dimension;
        radius = 1;
        tmin = 0.0;
        tmax = 1.0;
        V = (tmax - tmin) * volSphere(dimension - 1, radius);
    }

    TPoint* getRandomPoint() {
        double t = (tmin + (tmax - tmin) * rnd());
        vector<double> x = rndSpherical(this->dimension - 1, this->radius);
        return new Tp_EinsteinStaticUniverseD(this->dimension, t, x);
    }

    // spatial separation (length of the great circle segment)
    double arcDist(TPoint* _a, TPoint* _b) {
        Tp_EinsteinStaticUniverseD *a = (Tp_EinsteinStaticUniverseD*) _a;
        Tp_EinsteinStaticUniverseD *b = (Tp_EinsteinStaticUniverseD*) _b;
        if (a->dimension != b->dimension) throw "vectors not of same dimension";
        double sigma = 0.0;
        for (unsigned i = 0; i < a->dimension; i++) {
            sigma += a->x[i] * b->x[i];
        }
        return this->radius * acos(sigma);
    }

// causal relations
    bool prec(TPoint* _a, TPoint* _b) {
        Tp_EinsteinStaticUniverseD *a = (Tp_EinsteinStaticUniverseD*) _a;
        Tp_EinsteinStaticUniverseD *b = (Tp_EinsteinStaticUniverseD*) _b;
        return (b->t - a->t > arcDist(a, b));
    }

//	This gives the (modulus of the) geodesic distance between two events.
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        throw "not yet defined";
    }
    ;
};

