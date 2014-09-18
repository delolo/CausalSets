#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_Cylinder2 : public TPoint {
public:
	double x;
	Tp_Cylinder2(double _t, double _x) {
		t=_t;
		x=_x;
	}
	Tp_Cylinder2() {
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

// A geometry class for sprinkling onto a two dimensional flat Cylinder with circumference L and height T.
// Either into the full rectangle with width L and height T or into a causal interval between points
// (0,0) and (T,0).
class TGeometry_Cylinder2d : public TGeometry {
public:
	double L,T;
	bool causalInterval;
	Tp_Cylinder2 start,end;
	TGeometry_Cylinder2d(double _T, double _L) {
		T=_T;
		L=_L;
		causalInterval=false;
		V=T*L;
	}
	TGeometry_Cylinder2d(double _T, double _L, bool _causalInterval) {
		T=_T;
		L=_L;
		causalInterval=_causalInterval;
		if (causalInterval) {
			start=Tp_Cylinder2(0,0);
			end=Tp_Cylinder2(T,0);
			if(T>L) {
				V=0.5*T*T+(T-L)*L;
			} else {
				V=0.5*T*T;
			}
		} else {
			V=T*L;
		}
	}
	TPoint* getRandomPoint() {
		if (causalInterval) {
			Tp_Cylinder2 *candidate;
			bool inside=false;
			while (!inside) {
				candidate=new Tp_Cylinder2(rnd()*T,rnd()*L);
				if (prec(&start,candidate) && prec(candidate,&end)) {
					inside=true;
				} else {
					delete candidate;
				}
			}
			return candidate;
		} else {
			return new Tp_Cylinder2(rnd()*T,rnd()*L);
		}
	}
	double minDistance(Tp_Cylinder2 *a, Tp_Cylinder2 *b) {
		return min(abs(a->x-b->x),L-abs(a->x-b->x));
	}
	bool prec(TPoint *_a, TPoint *_b) {
		Tp_Cylinder2 *a=(Tp_Cylinder2*) _a;
		Tp_Cylinder2 *b=(Tp_Cylinder2*) _b;
		return ((b->t-a->t)>minDistance(a,b));
	}
	int countCausal(TPoint *_a, TPoint *_b) {
		Tp_Cylinder2 *a=(Tp_Cylinder2*) _a;
		Tp_Cylinder2 *b=(Tp_Cylinder2*) _b;
		if ((b->t-a->t)>minDistance(a,b)) {
			Tp_Cylinder2 *c=new Tp_Cylinder2(a->t+L/2,a->x+L/2);
			while (c->x>L) {
				c->x=c->x-L;
			}
			int Result=1+countCausal(c,b);
			delete c;
			return Result;
		} else {
			return 0;
		}
	}
};
