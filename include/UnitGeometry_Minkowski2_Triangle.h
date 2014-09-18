#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_Minkowski2_Rectangle : public TPoint {
public:
	double x;
	Tp_Minkowski2_Rectangle(double _t, double _x) {
		t=_t;
		x=_x;
	}
	Tp_Minkowski2_Rectangle() {
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

class TGeometry_Minkowski2_Triangle : public TGeometry {
public:
	double scale;
	TGeometry_Minkowski2_Triangle() {
		scale=1;
		V=scale*scale/2.;
	}
	TGeometry_Minkowski2_Triangle(double _scale) {
		scale=_scale;
		V=scale*scale/2.;
	}
	TPoint* getRandomPoint() {
		double u,v;
		bool valid=false;
		do {
			u=rnd()*scale;
			v=rnd()*scale;
			if (v<1-u) {
				valid=true;
			}
		} while (! valid);
		return new Tp_Minkowski2_Rectangle((u+v)/sqrt((double)2),(u-v)/sqrt((double)2));
	}
	bool prec(TPoint *_a, TPoint *_b) {
	    Tp_Minkowski2_Rectangle *a=(Tp_Minkowski2_Rectangle*) _a;
	    Tp_Minkowski2_Rectangle *b=(Tp_Minkowski2_Rectangle*) _b;
		return ((b->t-a->t)>abs(b->x-a->x));
	}
};
