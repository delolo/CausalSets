#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"
#include "UnitGeometry_MinkowskiD_Rectangle.h"


/*  A geometry class corresponding to a d-hypercube in spacetime with metric
 *
 *          ds^2 = (t/alpha)^2 (-dt^2 + dx^2).
 *
 *  The trace of the extrinsic curvature on a t=const. surface is given by
 *
 *          K = (d-1) * alpha / t^2
 *
 */
class TGeometry_CFTimeSquaredD_Rectangle: public TGeometry {
private:
    int dim;
    double alpha;
    double tmin, tmax;
    std::vector<double> xmax;
public:
    // constructor with specified side lengths:
    TGeometry_CFTimeSquaredD_Rectangle(int dimension, vector<double> xmax, double tmin, double tmax) {
        this->dim = dimension;
        this->alpha = 1.0;
        this->tmin = tmin;
        this->tmax = tmax;
        this->xmax = xmax;
        this->V = alpha / (dim + 1.0)
                * (pow(tmax / alpha, dim + 1) - pow(tmin / alpha, dim + 1));
        for (int i = 0; i < this->dim - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }
    // constructor with specified time interval:
    TGeometry_CFTimeSquaredD_Rectangle(int dimension, double tmin, double tmax) {
        this->dim = dimension;
        this->alpha = 1.0;
        this->tmin = tmin;
        this->tmax = tmax;
        this->xmax.resize(dimension - 1, 1.0);
        this->V = alpha / (dim + 1.0)
                * (pow(tmax / alpha, dim + 1) - pow(tmin / alpha, dim + 1));
        for (int i = 0; i < this->dim - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }
// constructor with unit side lengths:
    TGeometry_CFTimeSquaredD_Rectangle(int dimension) {
        this->alpha = 1.0;
        this->dim = dimension;
        this->tmin = 0.0;
        this->tmax = 1.0;
        this->xmax.resize(dimension - 1, 1.0);
        this->V = alpha / (dim + 1.0)
                * (pow(tmax / alpha, dim + 1) - pow(tmin / alpha, dim + 1));
        for (int i = 0; i < dim - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }

    void setAlpha(double alpha) {
        this->alpha = alpha;
        this->V = alpha / (dim + 1.0)
                * (pow(tmax / alpha, dim + 1) - pow(tmin / alpha, dim + 1));
        for (int i = 0; i < dim - 1; i++) {
            this->V = this->V * xmax[i];
        }
    }

    bool prec(TPoint *_a, TPoint *_b) {
        Tp_MinkowskiDCartesian *a = (Tp_MinkowskiDCartesian*) _a;
        Tp_MinkowskiDCartesian *b = (Tp_MinkowskiDCartesian*) _b;
        if (b->t < a->t) return false;
        double xdist = 0.0;
        for (int i = 0; i < this->dim - 1; i++) {
            xdist += sqr(b->x[i] - a->x[i]);
        }
        return (sqr(b->t - a->t) > xdist);
    }

    double geoDistanceReal(TPoint *_a, TPoint *_b) {
        return 0.0;
    }

    bool inSubregion(TPoint* _a) {
        return false;
    }

private:
    TPoint* getRandomPoint() {
        double n = pow(this->tmax / alpha, this->dim + 1)
                 - pow(this->tmin / alpha, this->dim + 1);
        double t = pow(this->tmin / alpha, this->dim + 1) + n * rnd();
        t = alpha * pow(t, 1 / (this->dim + 1.0));
        std::vector<double> x;
        for (int i = 0; i < dim - 1; i++) {
            x.push_back(this->xmax[i] * rnd());
        };
        return new Tp_MinkowskiDCartesian(dim, t, x);
    }
};
