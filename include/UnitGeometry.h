#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Tools.h"


// Defining the fundamental point class

// The class Tp and its derived classes hold the coordinates of a causal set element on the manifold.
class TPoint {
public:
	// The classes TSprinkling and TConnection both require
	// that in all spacetimes that will be dealt with there is a
	// global time coordinate t such that a point x can only be to
	// the causal past of y if x.t<y.t.
	// So all actual point classes inherit the coordinate t from the fundamental class.
	double t;
	TPoint(double _t) {
		t=_t;
	}
	TPoint() {
		t=0;
	}
	// The methods write and swrite are not necessary for causet computations.
	// They simply make writing causet data into files a little more convenient.
	// They are declared virtual so that the corresponding methods in the derived classes will be called instead.
	virtual void write(ofstream *os) {
		(*os) << t << endl;
	}
	virtual void swrite(ofstream *os) {
		(*os) << t;
	}
	virtual string toString() {
	    return "no toString() defined";
	}
	// virtual destructor (to avoid Eclipse complaints)
	virtual ~TPoint() {}
};

// Defining the fundamental geometry class

// The class TGeometry and its derived classes hold information about the manifold that shall be sprinkled in.
// All methods here are only virtual dummies that will be overwritten by actual Geometry implementations.
// See UnitGeometry_Minkowski2d or UnitGeometry_deSitter2d for examples.
class TGeometry {
public:
	// For some applications the volume of the sprinkling region needs to be known.
	// For example the d'Alembertian needs to know the Planck length of a sprinkling which is sqrt(V/N).
	double V;
	int SUBREGION;
	TGeometry() { V = 0.0; SUBREGION = 0; } // setting fields to default to avoid Eclipse warning.
	// getRandomPoint() will give back a pointer to a random sprinkling point.
	// Actual Geometry-implementations need to do a sprinkling map from the hypercube [0,1]^d to the manifold region.
	virtual TPoint* getRandomPoint() {
		return new TPoint();
	}
	// isCausal() decides whether a is to the causal past of b.
	virtual bool prec(TPoint *a, TPoint *b) {
		return false;
	}
	// geoDistanceReal() gives (the modulus of) the geodesic distance between a and b.
	// isCausal() decides whether a is to the causal past of b.
	virtual double geoDistanceReal(TPoint *a, TPoint *b) {
	    return 0;
	}
	// countCausal() gives the number of equivalence classes of causal trajectories from a to b.
	// This method only needs to be altered on geometries with nontrivial topology. Everywhere else there will
	// only be one equivalence class and thus this virtual implementation can be kept.
	// For a nontrivial example see UnitGeometry_Cylinder2d.
	virtual int countCausal(TPoint *a, TPoint *b) {
		return (int) prec(a,b);
	}
	// inSubregion() asserts whether a given element lies in a subregion of the geometry
	// as defined in the particular geometry class. For a nontrivial example see
	// UnitGeometry_Minkowski2_Rectangle.
	virtual bool inSubregion(TPoint* _a) {return false;}
	virtual ~TGeometry() {}
};
