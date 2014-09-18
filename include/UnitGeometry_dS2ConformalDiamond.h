#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_dS2Conformal : public TPoint {
public:
	double x;
	Tp_dS2Conformal(double _t, double _x) {
		t=_t;
		x=_x;
	}
	Tp_dS2Conformal() {
		t=0;
		x=0;
	}
	void write(ofstream *os) {
		(*os) << t << "\t" << x << endl;
	}
};

// A geometry class for sprinkling into a causal interval in deSitter 2d space in 
// conformally flat coordinates for times between (t0,0) and (t1,0)
class TGeometry_dS2ConformalDiamond : public TGeometry {
public:
	double alpha,t0,t1;
	TGeometry_dS2ConformalDiamond(double _alpha, double _t0, double _t1) {
		alpha=_alpha;
		t0=_t0;
		t1=_t1;
		V=2*sqr(alpha)*log(sqr(t0+t1)/(4*t0*t1));
	}
	TPoint* getRandomPoint() {
		double t,y;
		bool valid;
		// Points will be picked from a rectangle that contains the causal interval until a point
		// inside of the interval is picked.
		do {
			t=1/(1/t0-rnd()*(1/t0-1/t1));
			y=(rnd()-0.5)*(t1-t0);
			valid=true;
			if (t>(t0+t1)/2) {
				if (abs(y)>t1-t) {
					valid=false;
				}
			} else {
				if (abs(y)>t-t0) {
					valid=false;
				}
			}
		} while (! valid);
		return new Tp_dS2Conformal(t,y);
	}
	bool prec(TPoint *_a, TPoint *_b) {
		Tp_dS2Conformal *a=(Tp_dS2Conformal*) _a;
		Tp_dS2Conformal *b=(Tp_dS2Conformal*) _b;
		return (b->t-a->t>abs(b->x-a->x));
	}
//	This gives the (modulus of the) geodesic distance between two events. Since this is imaginary
//	for timelike separated events, the formula takes into account whether the events are causally
//	connected. CHECK THAT Z IS NEVER SMALLER THAN -1!!!
	double geoDistanceReal(TPoint *_a, TPoint *_b) {
		Tp_dS2Conformal *a=(Tp_dS2Conformal*) _a;
		Tp_dS2Conformal *b=(Tp_dS2Conformal*) _b;
		double z = 1.0+(sqr(b->t-a->t)-sqr(b->x-a->x))/(2.0*(a->t)*(b->t));
		if (z>1.0) {
			return alpha*log(z+sqrt(sqr(z)-1.0));
		} else {
			return alpha*acos(z);
		}
	}
};
