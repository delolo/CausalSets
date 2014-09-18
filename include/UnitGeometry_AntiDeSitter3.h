/*************************************************************************
 *
 *  The UnitGeometry class file for a finite region of 3-dimensional AdS.
 *  The coordinates (t,theta,phi) are for the metric
 *
 *  ds^2 = alpha^2 sec(theta)^2 (-dt^2 + dtheta^2 + sin(theta)^2 dphi^2),
 *
 *	where 	t 		is in (-infinity,infinity) for the universal cover,
 *			theta 	is in [0,PI/2) and
 *			phi		is in [0,2PI).
 *
 *	The sprinkling region correspond to t: 		[tmin,tmax]
 *										theta:	[0,thetamax]
 *										phi:	unrestricted.
 *
 *************************************************************************/

#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Tools.h"

class Tp_AntiDeSitter3: public TPoint {
public:
	double theta, phi;
	Tp_AntiDeSitter3(double _t, double _theta, double _phi) {
		t = _t;
		theta = _theta;
		phi = _phi;
	}
	Tp_AntiDeSitter3() {
		t = 0;
		theta = 0;
		phi = 0;
	}
	void write(ofstream *os) {
		(*os) << t << "\t" << theta << "\t" << phi << endl;
	}
};

//  A geometry class for sprinkling into a region of
//  the universal cover of ads2 in conformally flat coordinates
//  for times between 0 and 1
class TGeometry_AntiDeSitter3: public TGeometry {
public:
	double thetamax, tmin, tmax, alpha;
	// constructor with ads radius and time interval [0,1]
	TGeometry_AntiDeSitter3(double _thetamax, double _tmin, double _tmax) {
		alpha = 1.0;
		thetamax = _thetamax;
		tmin = _tmin;
		tmax = _tmax;
		V = alpha * alpha * PI * (pow(cos(thetamax), -2.0) - 1.0)
				* (tmax - tmin);
	}
// constructor without ads radius (set to =1) and time interval [0,1]
	TGeometry_AntiDeSitter3(double _thetamax) {
		alpha = 1.0;
		thetamax = _thetamax;
		tmin = 0.0;
		tmax = 1.0;
		V = alpha * alpha * PI * (pow(cos(thetamax), -2.0) - 1.0)
				* (tmax - tmin);
	}

	TPoint* getRandomPoint() {
		//TODO: put in radius alpha
		double t = alpha * (tmin + (tmax - tmin) * rnd());
		double theta = atan(sqrt(rnd() * (pow(cos(thetamax), -2.0) - 1.0)));
		double phi = 2 * PI * rnd();
		return new Tp_AntiDeSitter3(t, theta, phi);
	}

	// spatial separation (length of the great circle segment)
	// between two points in ads3
	double arcDist(TPoint* _a, TPoint* _b) {
		Tp_AntiDeSitter3 *a = (Tp_AntiDeSitter3*) _a;
		Tp_AntiDeSitter3 *b = (Tp_AntiDeSitter3*) _b;
		double d = cos(a->theta) * cos(b->theta)
				+ sin(a->theta) * sin(b->theta) * cos(b->phi - a->phi);
		return acos(d);
	}

	// does a precede b?
	bool prec(TPoint *a, TPoint *b) {
		return (b->t - a->t > arcDist(a, b));
	}
//	This gives the (modulus of the) geodesic distance between two events. Since this is imaginary
//	for timelike separated events, the formula takes into account whether the events are causally
//	connected. CHECK THAT Z IS NEVER SMALLER THAN -1!!!
	double geoDistanceReal(TPoint *_a, TPoint *_b) {
		throw "geoDistanceReal not yet defined for ads2";
		return 0.0;
	}

    // does a precede b after projection onto the boundary?
    bool precOnBoundary(TPoint* _a, TPoint* _b) {
        Tp_AntiDeSitter3 *a = (Tp_AntiDeSitter3*) _a;
        Tp_AntiDeSitter3 *b = (Tp_AntiDeSitter3*) _b;
        a->theta = PI / 2.0; // projection onto equator
        b->theta = PI / 2.0; // projection onto equator
        if (b->t - a->t >= arcDist(a, b)) return true;
        else return false;
    }
};

