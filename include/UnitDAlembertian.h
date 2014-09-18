#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
//#include <vector>
#include "UnitConnection.h"

using namespace std;

// TDAlembertian can be used to compute the two dimensional discrete d'Alembertian of a causal set.
// The constructor needs a reference to an instance of TConnectionWithIntervals that stores the causal matrix,
// the nonlocality scale lk and the discreteness scale l.

class TDAlembertian {
public:
	TConnectionWithIntervals *Connection;
	double lk,l,epsilon;
	TDAlembertian(TConnectionWithIntervals *_Connection, double _lk, double _l){
		Connection=_Connection;
		lk=_lk;
		l=_l;
		epsilon=pow(l/lk,2);
	}
	double f(int n) {
		double u=epsilon/(1-epsilon);
		return pow(1-epsilon,n)*(1+u*(-2.*n+u*n*(n-1)/2.));
	}
	double f_local(int n) {
		switch(n) {
			case 0: return 1;
			case 1: return -2;
			case 2: return 1;
			default: return 0;
		}
	}
	// This computes the nonlocal d'Alembertian evaluated at pos for a function phi
	// given by the vector such that phi_i = phi[i]
	double compute(vector<double> *phi, int pos) {
		double result=0;
		int a;
		for (a=0;a<pos;a++) {
			if (Connection->causal(a,pos)) {
				result+=f(Connection->cardinality(a,pos))*(*phi)[a];
			}
		}
		result=result*epsilon;
		result-=0.5*(*phi)[pos];
		result=4./pow(lk,2)*result;
		return result;
	}
	// This computes the nonlocal d'Alembertian for the constant function 1.
	double compute(int pos) {
		double result=0;
		int a;
		for (a=0;a<pos;a++) {
			if (Connection->causal(a,pos))
				result+=f(Connection->cardinality(a,pos));
		}
		result=result*epsilon;
		result-=0.5;
		result=4/pow(lk,2)*result;
		return result;
	}
	// This computes the local d'Alembertian for the constant function 1.
	double compute_local(int pos) {
		double result=0;
		int a;
		for (a=0;a<pos;a++) {
			if (Connection->causal(a,pos))
				result+=f_local(Connection->cardinality(a,pos));
		}
		result-=0.5;
		result=4/pow(l,2)*result;
		return result;
	}
};
