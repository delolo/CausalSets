#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_Minkowski2 : public TPoint {
public:
	double x;
	Tp_Minkowski2(double _t, double _x) {
		t=_t;
		x=_x;
	}
	Tp_Minkowski2() {
		t=0;
		x=0;
	}
	void write(ofstream *os) {
		(*os) << t << "\t" << x << endl;
	}
	void swrite(ofstream *os) {
		(*os) << t << "\t" << x;
	}
};

class TGeometry_Minkowski2 : public TGeometry {
public:
	double scale;
	TGeometry_Minkowski2() {
		scale=1;
		V=scale*scale;
	}
	TGeometry_Minkowski2(double _scale) {
		scale=_scale;
		V=scale*scale;
	}
	TPoint* getRandomPoint() {
		double u,v;
		u=rnd()*scale;
		v=rnd()*scale;
		return new Tp_Minkowski2((u+v)/sqrt((double)2),(u-v)/sqrt((double)2));
	}
	bool prec(TPoint *_a, TPoint *_b) {
		Tp_Minkowski2 *a=(Tp_Minkowski2*) _a;
		Tp_Minkowski2 *b=(Tp_Minkowski2*) _b;
		return ((b->t-a->t)>abs(b->x-a->x));
	}
	double geoDistanceReal(TPoint *_a, TPoint *_b) {
		Tp_Minkowski2 *a=(Tp_Minkowski2*) _a;
		Tp_Minkowski2 *b=(Tp_Minkowski2*) _b;
		if (abs(b->t-a->t)>abs(b->x-a->x)) {
			return sqrt(sqr(b->t-a->t)-sqr(b->x-a->x));
		} else {
			return sqrt(-sqr(b->t-a->t)+sqr(b->x-a->x));
		}
	}
};
