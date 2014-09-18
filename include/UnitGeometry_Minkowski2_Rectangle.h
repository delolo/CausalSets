#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_Minkowski2_Rectangle: public TPoint {
public:
    double x;
    Tp_Minkowski2_Rectangle(double _t, double _x) {
        t = _t;
        x = _x;
    }
    Tp_Minkowski2_Rectangle() {
        t = 0;
        x = 0;
    }
    void write(ofstream *os) {
        (*os) << t << "\t" << x << endl;
    }
    void swrite(ofstream *os) {
        (*os) << t << "\t" << x;
    }
//    string toString() {
//        return to_string(t) + "\t" + to_string(x);
//    }
};

class TGeometry_Minkowski2_Rectangle: public TGeometry {
public:
    double width;
    double height;
    /*  Integer specifying the type of subregion:
     *  0 = NO SUBREGION SET
     *  1 = RECTANGLE {x < x0}
     *  2 = RINDLER-LIKE PARTITION {(x + sub[1])^2-t^2 < sub[0]^2} */
    std::vector<double> sub;
    TGeometry_Minkowski2_Rectangle() {
        width = 1.0;
        height = 1.0;
        SUBREGION = 0;
        sub.resize(0);
        V = width * height;
    }
    TGeometry_Minkowski2_Rectangle(int _subregion, std::vector<double> _parameters) {
         width = 1.0;
         height = 1.0;
         V = width * height;
         SUBREGION = _subregion;
         sub = _parameters;
     }
    TGeometry_Minkowski2_Rectangle(double _height, double _width) {
        height = _height;
        width = _width;
        SUBREGION = 0;
        sub.resize(0);
        V = width * height;
    }
    TGeometry_Minkowski2_Rectangle(double _height, double _width, int _subregion, std::vector<double> _parameters) {
        height = _height;
        width = _width;
        SUBREGION = _subregion;
        sub = _parameters;
        V = width * height;
    }
    TPoint* getRandomPoint() {
        return new Tp_Minkowski2_Rectangle(height * rnd(), width * rnd());
    }
    bool prec(TPoint *_a, TPoint *_b) {
        Tp_Minkowski2_Rectangle *a = (Tp_Minkowski2_Rectangle*) _a;
        Tp_Minkowski2_Rectangle *b = (Tp_Minkowski2_Rectangle*) _b;
        return ((b->t - a->t) > abs(b->x - a->x));
    }
    bool inSubregion(TPoint* _a) {
        if(SUBREGION == 0) throw "NEED TO SPECIFY SUBREGION";
        if(SUBREGION == 1) return inSubregionRectangle(_a);
        if(SUBREGION == 2) return inSubregionRindler(_a);
        else return false;
    }
    bool inSubregionRectangle(TPoint* _a) {
        Tp_Minkowski2_Rectangle* a = (Tp_Minkowski2_Rectangle*) _a;
        if (a->x < this->sub[0]) return true;
        else return false;
    }
    bool inSubregionRindler(TPoint* _a) {
        Tp_Minkowski2_Rectangle* a = (Tp_Minkowski2_Rectangle*) _a;
        double t = a->t;
        double x = a->x;
        if (x < sqrt(t*t + sub[0]*sub[0]) -sub[1]) return true;
        else return false;
    }
};
