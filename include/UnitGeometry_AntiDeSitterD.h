/*************************************************************************
 *
 *  The UnitGeometry class file for a finite region of D-dimensional AdS.
 *
 *************************************************************************/

#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Tools.h"

class Tp_AntiDeSitterD: public TPoint {
public:
    int dimension;
    double t;
    std::vector<double> x;

    Tp_AntiDeSitterD(int _dimension, double _t, vector<double> _x) {
        if (_x.size() != dimension) throw "size of vector and dimension do not agree";
        dimension = _dimension;
        t = _t;
        x = _x;
    }
    Tp_AntiDeSitterD(double _dimension) {
        t = 0.0;
        dimension = _dimension;
        x.resize(dimension);
        fill(x.begin(), x.end(), 0.0);
        //TODO: not a valid point on the sphere
    }
    void write(ofstream *os) {
        (*os) << t << "\t";
        for (int i = 0; i < dimension; i++) {
            (*os) << x(i) << "\t";
        }
        (*os) << endl;
    }
};

//  A geometry class for sprinkling into a region of
//  the universal cover of ads2 in conformally flat coordinates
//  for times between 0 and 1
class TGeometry_AntiDeSitterD: public TGeometry {
public:
    int dimension;
    double thetamax, tmin, tmax, alpha;
    // constructor with ads radius and time interval [0,1]
    TGeometry_AntiDeSitterD(int _dimension, double _thetamax, double _tmin,
            double _tmax) {
        dimension = _dimension;
        alpha = 1.0;
        thetamax = _thetamax;
        tmin = _tmin;
        tmax = _tmax;
        //TODO V = alpha * alpha * PI * (pow(cos(thetamax), -2.0) - 1.0)
        //		* (tmax - tmin);
    }
// constructor without ads radius (set to =1) and time interval [0,1]
    TGeometry_AntiDeSitterD(int _dimension, double _thetamax) {
        dimension = _dimension;
        alpha = 1.0;
        thetamax = _thetamax;
        tmin = 0.0;
        tmax = 1.0;
        //TODO V = alpha * alpha * PI * (pow(cos(thetamax), -2.0) - 1.0)
        //		* (tmax - tmin);
    }

    TPoint* getRandomPoint() {
        //TODO: put in radius alpha
        double t = alpha * (tmin + (tmax - tmin) * rnd());
        vector<double> x = rndSpherical(dimension, alpha);
        if(x[0]<0) x[0] = - x[0];
        return new Tp_AntiDeSitterD(dimension, t, x);
    }

    // spatial separation (length of the great circle segment)
    // between two points in ads3
    double arcDist(TPoint* _a, TPoint* _b) {
        Tp_AntiDeSitterD *a = (Tp_AntiDeSitterD*) _a;
        Tp_AntiDeSitterD *b = (Tp_AntiDeSitterD*) _b;
        if (a->dimension!=b->dimension) throw "vectors not of same dimension";
        double sigma = 0.0;
        for(unsigned i = 0; i < a->dimension; i++){
          sigma += a->x[i] * b->x[i];
        }
        return acos(sigma);
    }

// causal relations
    bool prec(TPoint *a, TPoint *b) {
        return (b->t - a->t > arcDist(a, b));
    }

//	This gives the (modulus of the) geodesic distance between two events.
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        throw "not yet defined";
    }

    // does a precede b after projection onto the boundary?
    bool precOnBoundary(TPoint* _a, TPoint* _b) {
        throw "not yet defined";
     }
};

