#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_dS2Closed : public TPoint {
public:
	double x;
	Tp_dS2Closed(double _t, double _x) {
		t=_t;
		x=_x;
	}
	Tp_dS2Closed() {
		t=0;
		x=0;
	}
	void write(ofstream *os) {
		(*os) << t << "\t" << x << endl;
	}
};

// A geometry class for sprinkling into a slab in deSitter 2d space in 
// conformally flat coordinates for times between (t0,0) and (t1,0)
class TGeometry_dS2ClosedSlab : public TGeometry {
public:
	double alpha,t0,t1;
	TGeometry_dS2ClosedSlab(double _alpha, double _t0, double _t1) {
		alpha=_alpha;
		t0=_t0;
		t1=_t1;
		V=2*PI*alpha*alpha*(tan(t1)-tan(t0));
	}
    TGeometry_dS2ClosedSlab(double _alpha, double _t1) {
		t0=0.0;
        alpha=_alpha;
		t1=_t1;
		V=2*PI*alpha*alpha*(tan(t1)-tan(t0));
	}
	TPoint* getRandomPoint() {
		double t,y;
        y=rnd()*2*PI;
        t=atan(rnd()*(tan(t1)-tan(t0))+tan(t0));
		return new Tp_dS2Closed(t,y);
	}
	bool prec(TPoint *_a, TPoint *_b) {
		Tp_dS2Closed *a=(Tp_dS2Closed*) _a;
		Tp_dS2Closed *b=(Tp_dS2Closed*) _b;
		return (b->t-a->t>abs(b->x-a->x));
	}
    //	This gives the (modulus of the) geodesic distance between two events. Since this is imaginary
    //	for timelike separated events, the formula takes into account whether the events are causally
    //	connected.
	double geoDistanceReal(TPoint *_a, TPoint *_b) {
		Tp_dS2Closed *a=(Tp_dS2Closed*) _a;
		Tp_dS2Closed *b=(Tp_dS2Closed*) _b;
		double z = 1.0+(sqr(b->t-a->t)-sqr(b->x-a->x))/(2.0*(a->t)*(b->t));
		if (z>1.0) {
			return alpha*log(z+sqrt(sqr(z)-1.0));
		} else {
			return alpha*acos(z);
		}
	}
};
