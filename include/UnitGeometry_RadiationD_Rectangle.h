#pragma once
#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_RadiationD_Cartesian: public TPoint {
public:
    int dimension;
    std::vector<double> x;
    Tp_RadiationD_Cartesian(int dimension, double t, std::vector<double> x) {
        if (x.size() != dimension - 1) throw "size of vector and spatial dimension do not agree";
        this->dimension = dimension;
        this->t = t;
        this->x = x;
    }
    Tp_RadiationD_Cartesian(int dimension) {
        this->dimension = dimension;
        this->t = 0.0;
        fill(this->x.begin(), this->x.end(), 0.0);
    }
    void write(ofstream *os) {
        (*os) << t << "\t";
        for (int i = 0; i < dimension - 1; i++) {
            (*os) << x[i] << "\t";
        }
        (*os) << endl;
    }
    void swrite(ofstream *os) {
        (*os) << t << "\t";
        for (int i = 0; i < dimension - 1; i++) {
            (*os) << x[i] << "\t";
        }
    }
};

class TGeometry_RadiationD_Rectangle: public TGeometry {
public:
    int dimension;
    double tmin, tmax;
    std::vector<double> xmax;
    double alpha;
    // constructor with specified side lengths:
    TGeometry_RadiationD_Rectangle(int dimension, vector<double> xmax, double tmin, double tmax) {
        this->dimension = dimension;
        this->alpha = 1.0;
        this->tmin = tmin;
        this->tmax = tmax;
        this->xmax = xmax;
        this->V = 2 / (dimension + 1.0)
                * (pow(alpha * sqrt(tmax), dimension + 1)
                 - pow(alpha * sqrt(tmin), dimension + 1));
        for (int i = 0; i < dimension - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }
    // constructor with specified height:
    TGeometry_RadiationD_Rectangle(int dimension, double tmax) {
        this->dimension = dimension;
        this->alpha = 1.0;
        this->tmin = 0.0;
        this->tmax = tmax;
        this->xmax.resize(dimension - 1, 1.0);
        this->V = 2 / (dimension + 1.0)
                * (pow(alpha * sqrt(tmax), dimension + 1)
                 - pow(alpha * sqrt(tmin), dimension + 1));
        for (int i = 0; i < dimension - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }
// constructor with unit side lengths:
    TGeometry_RadiationD_Rectangle(int dimension) {
        this->dimension = dimension;
        this->alpha = 1.0;
        this->tmin = 0.0;
        this->tmax = 1.0;
        this->xmax.resize(dimension - 1, 1.0);
        this->V = 2.0 / (dimension + 1.0);
    }

    bool prec(TPoint *_a, TPoint *_b) {
        Tp_RadiationD_Cartesian *a = (Tp_RadiationD_Cartesian*) _a;
        Tp_RadiationD_Cartesian *b = (Tp_RadiationD_Cartesian*) _b;
        if (b->t < a->t) return false;
        double xdist = 0.0;
        for (int i = 0; i < this->dimension - 1; i++) {
            xdist += sqr(b->x[i] - a->x[i]);
        }
        return (sqr(2.0 * (sqrt(b->t) - sqrt(a->t)) / alpha) > xdist);
    }
    //TODO: not checked
    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        return 0.0;
    }
    //  TODO: make a virtual method in Geometry. don't know how to make arguments
    //  variable though.
    bool inSubregion(TPoint* _a) {
        return false;
    }

private:
    TPoint* getRandomPoint() {
        double delta = pow(this->tmax, (dimension + 1.0) / 2.0)
                     - pow(this->tmin, (dimension + 1.0) / 2.0);
        double t = 0.5 * pow(this->tmin, (dimension + 1.0) / 2.0) * rnd();
        t = pow(t, 2.0 / (dimension + 1.0));
        std::vector<double> x;
        for (int i = 0; i < dimension - 1; i++) {
            x.push_back(this->xmax[i] * rnd());
        };
        return new Tp_MinkowskiDCartesian(dimension, t, x);
    }
};
