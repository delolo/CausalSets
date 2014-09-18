#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_AntiDeSitter2: public TPoint {
public:
    double x;
    Tp_AntiDeSitter2(double _t, double _x) {
        t = _t;
        x = _x;
    }
    Tp_AntiDeSitter2() {
        t = 0;
        x = 0;
    }
    void write(ofstream *os) {
        (*os) << t << "\t" << x << endl;
    }
};

//  A geometry class for sprinkling into a finite slab of
//  the universal cover of ads2 in conformally flat coordinates
//  for times between 0 and 1
class TGeometry_AntiDeSitter2: public TGeometry {
public:
    double r1, t0, t1;
    // constructor with ads radius and time interval [0,1]
    TGeometry_AntiDeSitter2(double _r1, double _t0, double _t1) {
        r1 = _r1;
        t0 = _t0;
        t1 = _t1;
        V = 2.0 * tan(r1) * (t1 - t0);
    }
// constructor without ads radius (set to =1) and time interval [0,1]
    TGeometry_AntiDeSitter2(double _r1) {
        r1 = _r1;
        t0 = 0.0;
        t1 = 1.0;
        V = 2.0 * tan(r1);
    }

    TPoint* getRandomPoint() {
        return new Tp_AntiDeSitter2(t0 + (t1 - t0) * rnd(),
                atan(2 * tan(r1) * rnd() - tan(r1)));
    }
    // ads2 is conformal to minkowski2 so causal relations are easy
    bool prec(TPoint *_a, TPoint *_b) {
        Tp_AntiDeSitter2 *a = (Tp_AntiDeSitter2*) _a;
        Tp_AntiDeSitter2 *b = (Tp_AntiDeSitter2*) _b;
        return (b->t - a->t > abs(b->x - a->x));
    }
    // This gives the (modulus of the) geodesic distance between two events. Since this is imaginary
    // for timelike separated events, the formula takes into account whether the events are causally
    // connected. CHECK THAT Z IS NEVER SMALLER THAN -1!!!
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        throw "geoDistanceReal not yet defined for ads2";
        return 0.0;
    }
};

