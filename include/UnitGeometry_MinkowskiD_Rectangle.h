#pragma once
#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_MinkowskiDCartesian: public TPoint {
public:
    int dimension;
    std::vector<double> x;
    Tp_MinkowskiDCartesian(int dimension, double t, std::vector<double> x) {
        if (x.size() != dimension - 1) throw "size of vector and spatial dimension do not agree";
        this->dimension = dimension;
        this->t = t;
        this->x = x;
    }
    Tp_MinkowskiDCartesian(int dimension) {
        this->dimension = dimension;
        this->t = 0.0;
        fill(this->x.begin(), this->x.end(), 0.0);
    }
    void write(ofstream *os) {
        (*os) << this->toString() << endl;
    }
    void swrite(ofstream *os) {
        (*os) << this->toString();
    }
    string toString() {
        std::ostringstream ss;
        ss << std::fixed;
        ss << t;
        for (int i = 0; i < dimension - 1; i++) {
            ss << "\t" << x[i];
        }
        return ss.str();
    }
};

class TGeometry_MinkowskiD_Rectangle: public TGeometry {
public:
    int dimension;
    double tmax;
    std::vector<double> xmax;
    // constructor with specified side lengths:
    TGeometry_MinkowskiD_Rectangle(int dimension, vector<double> xmax, double tmax) {
        this->dimension = dimension;
        this->tmax = tmax;
        this->xmax = xmax;
        this->V = tmax;
        for (int i = 0; i < dimension - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }
    // constructor with specified height:
    TGeometry_MinkowskiD_Rectangle(int dimension, double tmax) {
        this->dimension = dimension;
        this->tmax = tmax;
        this->xmax.resize(dimension - 1, 1.0);
        this->V = tmax;
        for (int i = 0; i < dimension - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }
// constructor with unit side lengths:
    TGeometry_MinkowskiD_Rectangle(int dimension) {
        this->dimension = dimension;
        this->tmax = 1.0;
        this->xmax.resize(dimension - 1, 1.0);
        this->V = 1.0;
    }

    bool prec(TPoint *_a, TPoint *_b) {
        Tp_MinkowskiDCartesian *a = (Tp_MinkowskiDCartesian*) _a;
        Tp_MinkowskiDCartesian *b = (Tp_MinkowskiDCartesian*) _b;
        if (b->t < a->t) return false;
        double xdist = 0.0;
        for (int i = 0; i < this->dimension - 1; i++) {
            xdist += sqr(b->x[i] - a->x[i]);
        }
        return (sqr(b->t - a->t) > xdist);
    }
    //TODO: not checked
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        Tp_MinkowskiDCartesian *a = (Tp_MinkowskiDCartesian*) _a;
        Tp_MinkowskiDCartesian *b = (Tp_MinkowskiDCartesian*) _b;
        double xdist = 0.0;
        for (int i = 0; i < this->dimension - 1; i++) {
            xdist += sqr(b->x[i] - a->x[i]);
        }
        return sqrt(abs(sqr(b->t - a->t) - xdist));
    }
    //  TODO: make a virtual method in Geometry. don't know how to make arguments
    //  variable though.
    bool inSubregion(TPoint* _a) {
        return false;
    }

private:
    TPoint* getRandomPoint() {
        double t = this->tmax * rnd();
        std::vector<double> x;
        for (int i = 0; i < dimension - 1; i++) {
            x.push_back(this->xmax[i] * rnd());
        };
        return new Tp_MinkowskiDCartesian(dimension, t, x);
    }
};
