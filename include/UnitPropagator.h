#pragma once
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "Tools.h"
#include "UnitSprinkling.h"
//#include "UnitConnection.h"

using namespace std;

class TPropagator {
public:
	TSprinkling *Sprinkling;
	vector<vector<double> > storeDistance;
	vector<double> mainST;
	vector<double> causalST;
	vector<double> retardedST;
	double m, V;
	int size, utLength;
	TPropagator(TSprinkling *Sprinkling, double _m) {
		this->Sprinkling = Sprinkling;
		size = Sprinkling->size;
		utLength = size * (size + 1) / 2;
		m = _m;
		this->V = this->Sprinkling->geometry->V;
		storeDistance.resize(0);
		mainST.resize(0);
		causalST.resize(0);
		retardedST.resize(0);
	}

//	This creates an array of the (absolute value of the) geodesic distance between pairs of events.
//	The array is upper triangular but contains the information for all pairs as D(x,y)=D(y,x) and D(x,x)=0.
	void createDistanceFromSprinkling() {
		int i, j;
		storeDistance.resize(size);
		for (i = 0; i < size; i++) {
			storeDistance[i].resize(size - i - 1);
			for (j = 0; j < size - i - 1; j++) {
				storeDistance[i][j] = Sprinkling->geoDistanceReal(i, i + j + 1);
			}
		}
	}

//	This extracts the geodesic distance for any two points from the stored values.
	double geoDistance(int i, int j) {
		double zero = 0.0;
		if (i < j)
			return storeDistance[i][j - i - 1];
		else if (i > j)
			return storeDistance[j][i - j - 1];
		else
			return zero;
	}

//	This writes a matrix of (absolute values of) geodesic distances to an external file:
	void writeGeodesic(ofstream *os) {
		int i, j;
		cout << size << endl;
		for (i = 0; i < size; i++) {
			(*os) << geoDistance(i, 0);
			for (j = 1; j < size; j++) {
				(*os) << " " << geoDistance(i, j);
			}
			(*os) << endl;
		}
	}
//	This creates the unit upper triangular matrix (1+0.5aC) in stacked column-major format for use in LAPACK:

	void createMainST(TSprinkling *Sprinkling) {
		int i, j, pos = 1;
		double c = 0.5 * m * m * V / size;
		mainST.resize(utLength);
		mainST[0] = 1.0;
		for (i = 1; i < size; i++) {
			for (j = 0; j < i; j++) {
				mainST[pos] = c * Sprinkling->prec(j, i);
				pos++;
			}
			mainST[pos] = 1.0;
			pos++;
		}
	}

//	This creates the strictly upper triangular matrix C in stacked column-major format for
//	use in LAPACK:

	void createCausalST(TSprinkling *Sprinkling) {
		int i, j, pos = 1;
		causalST.resize(utLength);
		causalST[0] = 0.0;
		for (j = 1; j < size; j++) {
			for (i = 0; i < j; i++) {
				causalST[pos] = Sprinkling->prec(i, j);
				pos++;
			}
			causalST[pos] = 0.0;
			pos++;
		}
	}

//	This writes out a unit upper triangular in stacked format:
	void writeMainST() {

		for (int i = 0; i < utLength; i++) {
			cout << mainST[i] << endl;
		}
		cout << endl;
	}

	void writeRetardedST() {
		for (int i = 0; i < utLength; i++) {
			cout << retardedST[i] << endl;
		}
		cout << endl;
	}
//  calculate retarded propagator (taking O(N^3))
	double** getRetardedPropagator() {
		double c = -0.5;
		double b = -this->m * this->m
				/ (this->Sprinkling->size / this->Sprinkling->geometry->V);
		double** kRet = new double*[size];
		for (int i = 0; i < size; i++) {
			kRet[i] = new double[size];
		}
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++)
				kRet[i][j] = c * this->Sprinkling->prec(i, j);
		}

		if (size > 2) {
			for (int i = size - 3; i >= 0; i--) {
				for (int j = i + 2; j < size; j++) {
					for (int k = i + 1; k < j; k++) {
						kRet[i][j] += b * c
								* (double) this->Sprinkling->prec(i, j)
								* kRet[k][j];
					}
				}
			}
		}
		return kRet;
	}
};

//		duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
//		cout << "this took: " << duration << " seconds" << endl;

//		utcausal
//2. make coordinate arrays and possible bins
//3. define the continuum propagators in the respective geometry files.

// A simple implementation of the merge-sort algorithm is used to sort the points.
// It was chosen because it is easy to implement and scales like O(n log(n)).
