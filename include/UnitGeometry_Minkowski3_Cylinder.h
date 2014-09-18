#pragma once
#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_Minkowski3Polar: public TPoint {
public:
    double r, theta;
    Tp_Minkowski3Polar(double _t, double _r, double _theta) {
        t = _t;
        r = _r;
        theta = _theta;
    }
    Tp_Minkowski3Polar() {
        t = 0;
        r = 0;
        theta = 0;
    }
    void write(ofstream *os) {
        (*os) << t << "\t" << r << "\t" << theta << endl;
    }
    void swrite(ofstream *os) {
        (*os) << t << "\t" << r << "\t" << theta;
    }
};

class TGeometry_Minkowski3_Cylinder: public TGeometry {
public:
    double height, radius, subradius;
    TGeometry_Minkowski3_Cylinder() {
        this->radius = 1.0;
        this->subradius = radius;
        this->height = 1.0;
        this->V = 2 * PI * radius * height;
    }
    TGeometry_Minkowski3_Cylinder(double _subradius) {
        this->radius = 1.0;
        this->subradius = _subradius;
        this->height = 1.0;
        this->V = 2 * PI * radius * height;
    }
    TGeometry_Minkowski3_Cylinder(double _radius, double _height) {
        this->radius = _radius;
        this->subradius = radius;
        height = _height;
        V = 2 * PI * radius * height;
    }
    TGeometry_Minkowski3_Cylinder(double _radius, double _height, double _subradius) {
        this->radius = _radius;
        this->subradius = _subradius;
        height = _height;
        V = 2 * PI * radius * height;
    }
    bool prec(TPoint *_a, TPoint *_b) {
        Tp_Minkowski3Polar *a = (Tp_Minkowski3Polar*) _a;
        Tp_Minkowski3Polar *b = (Tp_Minkowski3Polar*) _b;
        if (b->t < a->t) return false;
        return (sqr(b->t - a->t)
                > sqr(b->r) + sqr(a->r)
                        - 2.0 * a->r * b->r * cos(a->theta - b->theta));
    }
    //TODO: not checked
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        Tp_Minkowski3Polar *a = (Tp_Minkowski3Polar*) _a;
        Tp_Minkowski3Polar *b = (Tp_Minkowski3Polar*) _b;
        double dsquared = sqr(b->t - a->t) - sqr(b->r) + sqr(a->r)
                + 2.0 * a->r * b->r * cos(a->theta - b->theta);
        return sqrt(abs(dsquared));
    }
    //  TODO: make a virtual method in Geometry. don't know how to make arguments
    //  variable though.
    bool inSubregion(TPoint* _a) {
        Tp_Minkowski3Polar* a = (Tp_Minkowski3Polar*) _a;
        if (a->r < this->subradius) return true;
        else return false;
    }
private:
    TPoint* getRandomPoint() {
        double theta = 2 * PI * rnd();
        double r = radius * sqrt(rnd());
        double t = height * rnd();
        return new Tp_Minkowski3Polar(t, r, theta);
    }
};
