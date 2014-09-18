#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_SpatialExponential2_Cartesian: public TPoint {
public:
    double x;
    Tp_SpatialExponential2_Cartesian(double _t, double _x) {
        t = _t;
        x = _x;
    }
    Tp_SpatialExponential2_Cartesian() {
        t = 0;
        x = 0;
    }
    void write(ofstream *os) {
        (*os) << t << "\t" << x << endl;
    }
    void swrite(ofstream *os) {
        (*os) << t << "\t" << x;
    }
};

class TGeometry_SpatialExponential2_Rectangle: public TGeometry {
public:
    double alpha;
    double width;
    double height;

    /*  Integer specifying the type of subregion:
     *  0 = NO SUBREGION SET
     *  1 = RECTANGLE {x < x0}
     *  2 = RINDLER-LIKE PARTITION {(x + sub[1])^2-t^2 < sub[0]^2} */
    std::vector<double> sub;

    TGeometry_SpatialExponential2_Rectangle() {
        alpha = 1.0;
        width = 1.0;
        height = 1.0;
        SUBREGION = 0;
        sub.resize(0);
        V = 0.5 * alpha * (exp(2.0 * width / alpha) - 1.0) * height;
    }
    TGeometry_SpatialExponential2_Rectangle(double _alpha) {
        alpha = _alpha;
        width = 1.0;
        height = 1.0;
        SUBREGION = 0;
        sub.resize(0);
        V = 0.5 * alpha * (exp(2.0 * width / alpha) - 1.0) * height;
    }
    TGeometry_SpatialExponential2_Rectangle(int _subregion, std::vector<double> _parameters) {
         alpha = 1.0;
         width = 1.0;
         height = 1.0;
         SUBREGION = _subregion;
         sub = _parameters;
         V = 0.5 * alpha * (exp(2.0 * width / alpha) - 1.0) * height;
     }
    TGeometry_SpatialExponential2_Rectangle(double _height, double _width, double _alpha) {
        alpha = _alpha;
        height = _height;
        width = _width;
        SUBREGION = 0;
        sub.resize(0);
        V = 0.5 * alpha * (exp(2.0 * width / alpha) - 1.0) * height;
    }
    TGeometry_SpatialExponential2_Rectangle(double _height, double _width, double _alpha, int _subregion, std::vector<double> _parameters) {
        alpha = _alpha;
        height = _height;
        width = _width;
        SUBREGION = _subregion;
        sub = _parameters;
        V = 0.5 * alpha * (exp(2.0 * width / alpha) - 1.0) * height;
    }
    TPoint* getRandomPoint() {
        double t = height * rnd();
        double x = 0.5 * alpha * log((exp(2.0 * width / alpha) * rnd() - 1.0) + 1.0);
        return new Tp_SpatialExponential2_Cartesian(t, x);
    }
    bool prec(TPoint *_a, TPoint *_b) {
        Tp_SpatialExponential2_Cartesian *a = (Tp_SpatialExponential2_Cartesian*) _a;
        Tp_SpatialExponential2_Cartesian *b = (Tp_SpatialExponential2_Cartesian*) _b;
        return ((b->t - a->t) > abs(b->x - a->x));
    }
    bool inSubregion(TPoint* _a) {
        if(SUBREGION == 0) throw "NEED TO SPECIFY SUBREGION";
        if(SUBREGION == 1) return inSubregionRectangle(_a);
        if(SUBREGION == 2) return inSubregionRindler(_a);
        else return false;
    }
    bool inSubregionRectangle(TPoint* _a) {
        Tp_SpatialExponential2_Cartesian* a = (Tp_SpatialExponential2_Cartesian*) _a;
        if (a->x < this->sub[0]) return true;
        else return false;
    }
    bool inSubregionRindler(TPoint* _a) {
        Tp_SpatialExponential2_Cartesian* a = (Tp_SpatialExponential2_Cartesian*) _a;
        double t = a->t;
        double x = a->x;
        if (x < sqrt(t*t+sub[0]*sub[0]) -sub[1]) return true;
        else return false;
    }
};
