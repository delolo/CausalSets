#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_dS2Conformal: public TPoint {
public:
    double x;
    Tp_dS2Conformal(double _t, double _x) {
        t = _t;
        x = _x;
    }
    Tp_dS2Conformal() {
        t = 0;
        x = 0;
    }
    void write(ofstream *os) {
        (*os) << t << "\t" << x << endl;
    }
};

//  A geometry class for sprinkling into a coordinate rectangle
//  [t0,t1] x [0, width] in the conformally flat chart of
//  two-dimensional de Sitter space
class TGeometry_dS2ConformalSlab: public TGeometry {
public:

    double alpha;
    double t0, t1;
    double width;

    TGeometry_dS2ConformalSlab(double _alpha, double _t0, double _t1, double _width) {
        alpha = _alpha;
        t0 = _t0;
        t1 = _t1;
        width = _width;
        V = alpha * alpha * width * (1 / t0 - 1 / t1);
    }
    TPoint* getRandomPoint() {
        double t, y;
        // Points are picked from a coordinate rectangle
        t = 1 / ((1 / t1 - 1 / t0) * rnd() + 1 / t0);
        y = width * rnd();
        return new Tp_dS2Conformal(t, y);
    }
    bool prec(TPoint *_a, TPoint *_b) {
        Tp_dS2Conformal *a = (Tp_dS2Conformal*) _a;
        Tp_dS2Conformal *b = (Tp_dS2Conformal*) _b;
        return (b->t - a->t > abs(b->x - a->x));
    }
//  This gives the (modulus of the) geodesic distance between two events. Since this is imaginary
//  for timelike separated events, the formula takes into account whether the events are causally
//  connected.
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        Tp_dS2Conformal *a = (Tp_dS2Conformal*) _a;
        Tp_dS2Conformal *b = (Tp_dS2Conformal*) _b;
        double z = 1.0
                + (sqr(b->t - a->t) - sqr(b->x - a->x))
                        / (2.0 * (a->t) * (b->t));
        if (z > 1.0) {
            return alpha * log(z + sqrt(sqr(z) - 1.0));
        }
        else {
            return alpha * acos(z);
        }
    }
};
