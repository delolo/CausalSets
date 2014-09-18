/*----------------------------------------------------------------------------
 *
 *  Written:       10/10/2013
 *  Last updated:  24/02/2014
 *
 *  TConnection is handling the adjacency matrices A_C and A_R of the causal
 *  set. Some applications make excessive use of the causal relations of
 *  almost all points. Thus computing them once at the start might save some
 *  time. In particular applications that need to know the links the class
 *  TConnection is vital. As all points are assumed to be sorted by a global
 *  time coordinate t such that b can only be to the future of a if b.t>a.t
 *
 *  This means the adjacency matrices are upper triangular. So to save storage
 *  space not the full rectangle is stored but only the upper triangle
 *  storeCausal has size N, storeCausal[i] has size N-i-1 and storeCausal[i][j]
 *  stores whether element i is to the causal past of i+j+1.
 *
 *  The function causal(i,j) can be treated like the full square causal matrix.
 *  The very same is true for links.
 *
 *---------------------------------------------------------------------------*/

#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "UnitSprinkling.h"
#include <sstream>
#include "Tools.h"

using namespace std;

class TConnection {
public:

    vector<vector<bool> > causalmatrix;
    vector<vector<bool> > linkmatrix;
    int size;
    int subsize;
    vector<bool> atomsInSubregion;

    // Basic Constructor:
    TConnection(TSprinkling *Sprinkling) {
        atomsInSubregion.resize(0);
        causalmatrix.resize(0);
        linkmatrix.resize(0);
        size = 0;
        subsize = 0;
        createAdjacencyMatrixFromSprinkling(Sprinkling);
    }

    // The link matrix requires a lot of effort to be computed.
    // Thus for applications where it is not required it can be skipped by
    // calling with needLinks=false.
    //
    TConnection(TSprinkling *Sprinkling, bool needLinks, bool needSubregion) {
        atomsInSubregion.resize(0);
        causalmatrix.resize(0);
        linkmatrix.resize(0);
        size = 0;
        subsize = 0;
        createAdjacencyMatrixFromSprinkling(Sprinkling);
        if (needSubregion) {
            findAtomsInSubregion(Sprinkling);
        }
        if (needLinks) {
            createLinkMatrix();
        }
    }

    /* BASIC CAUSAL RELATIONS METHODS */

    // gives the total number of (irreflexive) relations in the causet:
    int relations() {
        int out = 0;
        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                if (this->causal(i, j)) out++;
            }
        }
        return out;
    }

    // This creates the upper triangular causal matrix and stores it in storageCausal.
    // Cost is O(N^2).
    void createAdjacencyMatrixFromSprinkling(TSprinkling *Sprinkling) {
        int i, j;
        size = Sprinkling->size;
        causalmatrix.resize(size);
        for (i = 0; i < size; i++) {
            causalmatrix[i].resize(size - i - 1);
            for (j = 0; j < size - i - 1; j++) {
                causalmatrix[i][j] = Sprinkling->prec(i, i + j + 1);
            }
        }
    }

    /*  For applications it is much more convenient to have something that looks more
     *  like a full matrix. For this the function causal(i,j) can be used. causal(i,j)
     *  returns whether element i in the Sprinkling->vertex is to the causal past
     *  of element j. */
    bool causal(int i, int j) {
        if (j <= i)
        return false;
        return causalmatrix[i][j - i - 1];
    }

    // To decide whether two causally related elements are linked one must go over all
    // potential elements in between and rule them out. Cost is O(N^3).
    void createLinkMatrix() {
        int i, j, k;
        linkmatrix.resize(size);
        for (i = 0; i < size; i++) {
            linkmatrix[i].resize(size - i - 1);
            for (j = i + 1; j < size; j++) {
                if (causal(i, j)) {
                    linkmatrix[i][j - i - 1] = true;
                    k = i + 1;
                    while ((k < j) && linkmatrix[i][j - i - 1]) {
                        if (causal(i, k) && causal(k, j)) {
                            linkmatrix[i][j - i - 1] = false;
                        }
                        k++;
                    }
                }
            }
        }
    }

    /*  See causal(i,j). This is the same for links. */
    bool link(int i, int j) {
        if (j <= i) return false;
        if (linkmatrix.size() == 0) return false;
        return linkmatrix[i][j - i - 1];
    }

    /* FANCIER CAUSAL RELATIONS METHODS */

    int cardOfFutureSet(int i) {
        if (i == this->size - 1) return 0;
        int out = 0;
        for (int j = i + 1; j < this->size; j++) {
            if (causal(i,j)) out++;
        }
        return out;
    }

    int cardOfPastSet(int i) {
        if (i == 0) return 0;
        int out = 0;
        for (int j = 0; j < i; j++) {
            if (causal(j,i)) out++;
        }
        return out;
    }


   /* SUBREGION METHODS */

    // this finds the elements in the subregion specified in a Sprinkling->geometry
    void findAtomsInSubregion(TSprinkling* Sprinkling) {
        atomsInSubregion.resize(Sprinkling->size);
        for (int i = 0; i < Sprinkling->size; i++) {
            atomsInSubregion[i] = Sprinkling->inSubregion(i);
            if(inSubregion(i)) subsize++;
        }
    }

    bool inSubregion(int i) {
        if(atomsInSubregion.size()==0) throw "Need to initiate atomsInSubregion to use this";
        return atomsInSubregion[i];
    }

    /* IO STUFF */

    // adjacency matrix as a string
    string toString()
    {
        std::ostringstream str;
        str << "Causal Matrix:\n";
        int i, j;
        for (i = 0; i < size; i++) {
            str << (int) causal(i, 0);
            for (j = 1; j < size; j++) {
                str << " " << (int) (causal(i, j));
            }
            str << "\n";
        }
        return str.str();
    }

    // This writes the causal matrix to a file.
    // The files will usually be really big and not very handy.
    // They can be used however to import the matrix into
    // contemporary computer algebra systems to do
    // things like propagator computation there.
    void writeCausalMatrix(string filename) {
        ofstream os;
        os.open(filename.c_str());
        int i, j;
        for (i = 0; i < size; i++) {
            os << (int) causal(i, 0);
            for (j = 1; j < size; j++) {
                os << " " << (int) (causal(i, j));
            }
            os << endl;
        }
        cout << "Printed causal matrix to the file \"" << filename << "\"."
                << endl;
        os.close();
    }

    // This writes the causal matrix to a file.
    // The files will usually be really big and not very handy.
    // They can be used however to import the matrix into
    // contemporary computer algebra systems to do
    // things like propagator computation there.
    void writeCausalList(string filename) {
        ofstream os;
        os.open(filename.c_str());
        int i, j;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if(causal(i,j)) os << i << "\t" << j << "\n";
            }
        }
        cout << "Printed causal list to the file \"" << filename << "\"."
                << endl;
        os.close();
    }

    // Same for link matrix.
    void writeLinkMatrix(ofstream *os) {
        int i, j;
        for (i = 0; i < size; i++) {
            (*os) << (int) link(i, 0);
            for (j = 1; j < size; j++) {
                (*os) << " " << (int) link(i, j);
            }
            (*os) << endl;
        }
    }
};

// THomotopyMatrix is literally identical to TConnection without the link part.
// Except that the usual causal relation is replaced by the relation that counts
// the number of equivalence classes of causal paths. This can be used to compute
// for example the propagator on a cylinder.

class THomotopyMatrix {
public:
    vector<vector<int> > storeCountCausal;
    int size;
    THomotopyMatrix(TSprinkling *Sprinkling) {
        storeCountCausal.resize(0);
        size = 0;
        createCausalFromSprinkling(Sprinkling);
    }
    void createCausalFromSprinkling(TSprinkling *Sprinkling) {
        int i, j;
        size = Sprinkling->size;
        storeCountCausal.resize(size);
        for (i = 0; i < size; i++) {
            storeCountCausal[i].resize(size - i - 1);
            for (j = 0; j < size - i - 1; j++) {
                storeCountCausal[i][j] = Sprinkling->countCausal(i, i + j + 1);
            }
        }
    }
    int countCausal(int i, int j) {
        if (j <= i)
        return false;
        return storeCountCausal[i][j - i - 1];
    }
    void writecountCausal(string filename) {
        ofstream os;
        os.open(filename.c_str());
        int i, j;
        for (i = 0; i < size; i++) {
            os << countCausal(i, 0);
            for (j = 1; j < size; j++) {
                os << " " << countCausal(i, j);
            }
            os << endl;
        }
        os.close();
    }
};
