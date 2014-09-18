#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "UnitConnection.h"

using namespace std;

// THomotopyAlgorithm implements the algorithm introduced in Chapter 4 of the thesis in both the original
// not so well working version and the improved version.
// It will try to determine which homotopy zones relative to the point given by _origin in the constructor all other set
// elements are in. count[i] will hold the number of the homotopy zone of element i with respect to origin.
// First all elements to the future of origin are initialized with 1, all others with 0.
// Then the algorithm will at each step determine the set Omega_origin(n) by either calling
// getOmega: original version that looks for minimal elements in H_origin^+(n)
// getOmegaMax: modified version that looks for maximal elements in \mathcal{C} \setminus H_origin^+(n).
// Both getOmega methods will write a list of points they find into the vector<int> omega that is passed along.
// liftZone will increase count[i] for all i that are to the future of all elements in omega. Thus it determines the set
// H_origin^+(n+1).

class THomotopyAlgorithm {
	public:
	TConnectionWithIntervals *Connection;
	int origin;
	static const int omegaDepth=0;
	vector<int> count;
	THomotopyAlgorithm(int _origin, TConnectionWithIntervals *_Connection) {
		origin=_origin;
		Connection=_Connection;
		constructCount(true);
	}
	THomotopyAlgorithm(int _origin, TConnectionWithIntervals *_Connection, bool maximal) {
		origin=_origin;
		Connection=_Connection;
		constructCount(maximal);
	}
	void constructCount(bool maximal) {
		count.resize(Connection->size);
		int i;
		for(i=0;i<Connection->size;i++) {
			count[i]=(int)(Connection->causal(origin,i));
		}
		vector<int> omega;
		bool cont=true;
		int step=1;
		while (cont) {
			if (maximal) {
				getOmegaMax(&omega,step);
			} else {
				getOmega(&omega,step);
			}
			cont=liftZone(&omega,step);
			step++;
		}
	}
	void getOmega(vector<int> *omega, int step) {
		omega->resize(0);
		bool isInOmega;
		int depthCounter;
		int i,j;
		for (i=origin+1;i<Connection->size;i++) {
			if((count[i]==step) && (i!=origin)) {
				isInOmega=true;
				depthCounter=0;
				j=origin+1;
				while((j<i) && isInOmega) {
					if((count[j]==step) && Connection->causal(j,i)) {
						depthCounter++;
						if(depthCounter>omegaDepth) {
							isInOmega=false;
						}
					}
					j++;
				}
				if (isInOmega) {
					omega->resize(omega->size()+1);
					(*omega)[omega->size()-1]=i;
				}
			}
		}
	}
	void getOmegaMax(vector<int> *omega, int step) {
		omega->resize(0);
		bool isInOmega;
		int depthCounter;
		int i,j;
		for (i=0;i<Connection->size;i++) {
			if((count[i]==step-1)) {
				isInOmega=true;
				depthCounter=0;
				j=i+1;
				while((j<Connection->size) && isInOmega) {
					if((count[j]==step-1) && Connection->causal(i,j)) {
						depthCounter++;
						if(depthCounter>omegaDepth) {
							isInOmega=false;
						}
					}
					j++;
				}
				if (isInOmega) {
					omega->resize(omega->size()+1);
					(*omega)[omega->size()-1]=i;
				}
			}
		}
	}
	bool liftZone(vector<int> *omega, int step) {
		int i,j;
		bool isInNextZone;
		bool foundOne=false;
		for(j=origin+1;j<Connection->size;j++) {
			if(count[j]==step){
				i=0;
				isInNextZone=true;
				while((i<(int)omega->size()) && isInNextZone) {
					if (! Connection->causal((*omega)[i],j)) {
						isInNextZone=false;
					}
					i++;
				}
				if(isInNextZone) {
					count[j]=step+1;
					foundOne=true;
				}
			}
		}
		return foundOne;
	}
};