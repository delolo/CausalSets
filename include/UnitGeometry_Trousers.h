#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include"Tools.h"

class Tp_Trousers : public TPoint {
public:
	double x;
	Tp_Trousers(double _t, double _x) {
		t=_t;
		x=_x;
	}
	Tp_Trousers() {
		t=0;
		x=0;
	}
	void write(ofstream *os) {
		(*os) << t << "\t" << x << endl;
	}
};

class TGeometry_Trousers : public TGeometry {
public:
	double z,t0,t1;
	TGeometry_Trousers(double _z, double _t0, double _t1) {
		z=_z;
		t0=_t0;
		t1=_t1;
		V=t1;
	}
	TPoint* getRandomPoint() {
		return new Tp_Trousers(rnd()*t1,rnd());
	}
	double minDistance(double ax, double bx, double circ) {
		return min(abs(ax-bx),circ-abs(ax-bx));
	}
	bool prec(TPoint *_a, TPoint *_b) {
		return prec(_a,_b,false);
	}
	bool prec(TPoint *_a, TPoint *_b, bool verbose) {
		Tp_Trousers *a=(Tp_Trousers*) _a;
		Tp_Trousers *b=(Tp_Trousers*) _b;
		if ((a->t<t0)&&(b->t<t0)) {
			if ((a->x<z) && (b->x<z)) {
				if (b->t-a->t>minDistance(a->x,b->x,z)) {
					return true;
				} else {
					return false;
				}
			}
			if ((a->x>z) && (b->x>z)) {
				if (b->t-a->t>minDistance(a->x-z,b->x-z,1-z)) {
					return true;
				} else {
					return false;
				}
			}
			return false;
		}
		if ((a->t>t0)&&(b->t>t0)) {
			if (b->t-a->t>minDistance(a->x,b->x,1)) {
				return true;
			} else {
				return false;
			}
		}
		if ((b->t>t0)&&(a->t<t0)) {
			double al,ar,bl,br,width,ax;
			// determine which part of the cylinder merging area lies in the future of a
			if (a->x<z) {
				ax=a->x;
				width=z;
			} else {
				ax=a->x-z;
				width=1-z;
			}
			al=ax-(t0-a->t);
			ar=ax+(t0-a->t);
			while(al<0) {
				al+=width;
				ar+=width;
			}
			if (ar-al>=width) {
				al=0;
				ar=width;
			}
			if (verbose) {
				cout << al << " " << ar << endl;
			}
			// determine which part of the cylinder merging area lies in the past of b
			// first determine the past of b on the broad upper cylinder part
			bl=b->x-(b->t-t0);
			br=b->x+(b->t-t0);
			while(bl<0) {
				bl+=1;
				br+=1;
			}
			if(br-bl>=1) {
				bl=0;
				br=1;
			}
			// now map the waist to the approriate leg
			if (a->x<z) {
				if(bl<z) {
					if(br<z) {
					} else {
						if(br<1) {
							br=z;
						} else {
							br-=(1-z);
						}
					}
				} else {
					bl=0;
					if(br<1) {
						//br=0;
						return false;
					} else {
						if (br<1+z) {
							br-=1;
						} else {
							if(br>1+z) {
								br=z;
							}
						}
					}
				}
			} else {
				if(bl<z) {
					bl=0;
					if(br<z) {
						br=0;
					} else {
						if(br<1) {
							br-=z;
						} else {
							br=(1-z);
						}
					}
				} else {
					bl-=z;
					if(br<1) {
						br-=z;
					} else {
						if(br<1+z) {
							br=(1-z);
						} else {
							br-=(2*z);
						}
					}
				}
			}
			// now compare if both cylinder areas overlap
			if ((al<br) && (ar>bl)) {
				return true;
			}
			if ((al+width<br) && (ar+width>bl)) {
				return true;
			}
			if ((al<br+width) && (ar>bl+width)) {
				return true;
			}
		}
		return false;
	}
};
