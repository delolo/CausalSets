#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include"UnitGeometry.h"

using namespace std;

// TSprinkling is creating and handling the list of sprinkled points into a given geometry.
// A list of pointers to the chosen spacetime points is stored in the vector vertex.
// The whole concept of Koko is the following:
// The fundamental classes Tp and TGeometry are defined in UnitGeometry.
// They contain no information about a real spacetime whatsoever. But they have virtual placeholders
// for all vital things a real geometry must supply when implemented.
// Like this TSprinkling works independent of the actual geometry.
// A real geometry will then be defined as a derived class of TGeometry and like this
// compatibility with TSprinkling is ensured.

class TSprinkling {
public:
    vector<TPoint*> points;
    vector<TPoint*> subvertex;
    TGeometry *geometry;
    int size;
    // The basic constructor with expected number of elements specified.
    // cont_exp gives: the number of expected elements
    // _geometry: pointer to an instance of the desired geometry.
    TSprinkling(int count_exp, TGeometry *_geometry) {
        // The number of points for a given density is picked by a Poisson distribution.
        // The Poisson distribution is rather ugly on a computer. For large count_exp it is
        // approximated very good by a Gaussian with same mean and width.
        geometry = _geometry;
        int count = round(abs(rndGaussian(count_exp, sqrt((double) count_exp))));
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        addRandom(count);
        sort();
    }
    //  The basic constructor with density specified.
    TSprinkling(double rhoexp, TGeometry *_geometry) {
        // The number of points for a given density is picked by a Poisson distribution.
        // The Poisson distribution is rather ugly on a computer. For large count_exp it is
        // approximated very good by a Gaussian with same mean and width.
        geometry = _geometry;
        int countexp = round(rhoexp * geometry->V);
        int count = round(abs(rndGaussian(countexp, sqrt((double) countexp))));
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        addRandom(count);
        sort();
    }
    //  The basic constructor with rho = 1.0;
    TSprinkling(TGeometry *_geometry) {
        // The number of points for a given density is picked by a Poisson distribution.
        // The Poisson distribution is rather ugly on a computer. For large count_exp it is
        // approximated very good by a Gaussian with same mean and width.
        geometry = _geometry;
        int countexp = geometry->V;
        int count = round(abs(rndGaussian(countexp, sqrt((double) countexp))));
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        addRandom(count);
        sort();
    }

    TSprinkling(double rho_exp, bool randomCount, TGeometry *_geometry) {
        // TODO: Need to check if a uniform sprinkling with rho ~ Poiss(rho_exp)
        // gives a Poisson sprinkling.
        geometry = _geometry;
        int count = 0;
        if (randomCount) {
            int count_exp = round(rho_exp * geometry->V);
            count = round(abs(rndGaussian(count_exp, sqrt((double) count_exp))));
        }
        else {
            count = round(rho_exp * geometry->V);
        }
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        addRandom(count);
        sort();
    }
    // Sometimes one might want the number of elements to be fixed. This can be done here by calling with
    // randomCont=false.
    TSprinkling(int count_exp, bool randomCount, TGeometry *_geometry) {
        int count = 0;
        if (randomCount) {
            count = round(abs(rndGaussian(count_exp, sqrt((double) count_exp))));
        }
        else {
            count = count_exp;
        }
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        geometry = _geometry;
        addRandom(count);
        sort();
    }
    // If a given point must be part of the sprinkling it can be given by a.
    TSprinkling(int count_exp, TPoint *a, TGeometry *_geometry) {
        int count = round(abs(rndGaussian(count_exp, sqrt((double) count_exp))));
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        geometry = _geometry;
        add(a);
        addRandom(count - 1);
        sort();
    }
    // Combining the last two constructors.
    TSprinkling(int count_exp, bool randomCount, TPoint *a, TGeometry *_geometry) {
        int count;
        if (randomCount) {
            count = round(abs(rndGaussian(count_exp, sqrt((double) count_exp))));
        }
        else {
            count = count_exp;
        }
        points.resize(0);
        subvertex.resize(0);
        size = 0;
        geometry = _geometry;
        add(a);
        addRandom(count - 1);
        sort();
    }
    // The destructor cleans up by deleting all the sprinkling points.
    ~TSprinkling() {
        for (int i = 0; i < points.size(); i++) {
            delete points[i];
        }
        if (!subvertex.empty()) {
            for (int i = 0; i < subvertex.size(); i++) {
                delete subvertex[i];
            }
        }
        vector<TPoint*>().swap(points);
        vector<TPoint*>().swap(subvertex);
    }
    // count points will be added to the sprinkling by requesting them from the geometry. geometry will deal
    // with the choice of the points. See UnitGeometry for more information.
    void addRandom(int count) {
        int i;
        points.resize(size + count);
        for (i = 0; i < count; i++) {
            points[size + i] = geometry->getRandomPoint();
        }
        size += count;
    }
    void add(TPoint *point) {
        size++;
        points.resize(size);
        points[size - 1] = point;
    }
    void remove(int i) {
        std::vector<TPoint*>::iterator it;
        it = points.begin() + i;
        delete * it;
        points.erase(it);
        this->size--;
    }
    // The idea is that "above" TSprinkling no other code parts shall interact with the TGeometry class directly.
    // So the causal relations functions are simply handed downwards by the TSprinkling class.
    // Here prec can be either referred to directly with pointers to the points in question or with their indices in the
    // vertex list.
    bool prec(TPoint *a, TPoint *b) {
        return geometry->prec(a, b);
    }
    bool prec(int i, int j) {
        return prec(points[i], points[j]);
    }
    // countCausal will give the number of equivalence classes of causally accessible paths between a and b.
    // This is only relevant for manifolds with nontrivial topology.
    int countCausal(TPoint *a, TPoint *b) {
        return geometry->countCausal(a, b);
    }
    int countCausal(int i, int j) {
        return countCausal(points[i], points[j]);
    }

    // Here we retrieve the geodesic distance function from the geometry file:
    double geoDistanceReal(TPoint *a, TPoint *b) {
        return geometry->geoDistanceReal(a, b);
    }
    double geoDistanceReal(int i, int j) {
        return geoDistanceReal(points[i], points[j]);
    }
    bool inSubregion(int i) {
        return geometry->inSubregion(points[i]);
    }
    void getSubvertex() {
        for (int i = 0; i < this->size; i++) {
            if (inSubregion(i)) subvertex.push_back(points[i]);
        }
    }
    // To save some computation effort all points are sorted by their t coordinate
    // which is inherited from the fundamental Tp class. So only points b with b.t>a.t
    // can be to the future of a.
    void sort() {
        mergeSort(&points);
    }
    // A simple implementation of the merge-sort algorithm is used to sort the points.
    // It was chosen because it is easy to implement and scales like O(n log(n)).
    void mergeSort(vector<TPoint*> *list) {
        if (list->size() <= 1)
        return;
        if (list->size() == 2) {
            if ((*list)[0]->t > (*list)[1]->t) {
                TPoint* buffer = (*list)[0];
                (*list)[0] = (*list)[1];
                (*list)[1] = buffer;
            }
            return;
        }
        int a = list->size() / 2;
        int i;
        vector<TPoint*> partA(a);
        vector<TPoint*> partB(list->size() - a);
        for (i = 0; i < a; i++) {
            partA[i] = (*list)[i];
        }
        for (i = a; i < (int) list->size(); i++) {
            partB[i - a] = (*list)[i];
        }
        mergeSort(&partA);
        mergeSort(&partB);
        int posA = 0, posB = 0, pos = 0;
        while ((pos < (int) list->size()) && (posA < (int) partA.size())
                && (posB < (int) partB.size())) {
            if (partA[posA]->t < partB[posB]->t) {
                (*list)[pos] = partA[posA];
                posA++;
            }
            else {
                (*list)[pos] = partB[posB];
                posB++;
            }
            pos++;
        }
        for (i = posA; i < (int) partA.size(); i++) {
            (*list)[pos] = partA[i];
            pos++;
        }
        for (i = posB; i < (int) partB.size(); i++) {
            (*list)[pos] = partB[i];
            pos++;
        }
    }

    string toString() {
        std::ostringstream ss;
        for (int i = 0; i < this->size; i++) {
            ss << points[i]->toString() << "\n";
        }
        return ss.str();
    }

    string stringOfCausal() {
        string str;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (prec(i, j)) str += "1";
                else str += "0";
                str += " ";
            }
            str += "\n";
        }
        return str;
    }
    // Writes a list of all point coordinates to a file.
    void writeToFile(string filename) {
        ofstream os;
        os.open(filename.c_str());
        int i;
        for (i = 0; i < size; i++) {
            points[i]->write(&os);
        }
        cout << "Printed the coordinates of " << this->size
                << " elements to the file \"" << filename
                << "\"." << endl;
        os.close();
    }
    // Writes a list of all point coordinates to a file
    // with boolean inSubregion
    void writeToFileWithSubregion(string filename) {
        ofstream os;
        os.open(filename.c_str());
        int i;
        for (i = 0; i < size; i++) {
            os << inSubregion(i) << "\t";
            points[i]->write(&os);
        }
        cout << "Printed the coordinates of " << this->size
                << " elements to the file \"" << filename
                << "\"." << endl;
        os.close();
    }
};
