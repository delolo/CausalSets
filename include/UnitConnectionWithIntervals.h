/*----------------------------------------------------------------------------
 *
 *  Written:       10/10/2013
 *  Last updated:  24/02/2014
 *
 *  This is an extension of TConnection that contains methods for finding
 *  intervals, chains, and so on.
 *
 *  The methods for intervals scale at most O(N^3). The brute force method to
 *  find the number of card-n intervals is to square A_C and count how many
 *  entries are equal to n. But squaring A_C takes at least O(N^2.8) using the
 *  Strassen algorithm. Instead, the methods that count card-n order
 *  intervals have been optimized to give better scaling by only counting
 *  intervals of card-n exactly (and not those of larger cardinality). The
 *  methods for N1X, N1Y and D1 scale roughly like O(N^2.3).
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

class TConnectionWithIntervals {
public:

    vector<vector<bool> > causalmatrix;
    vector<vector<bool> > linkmatrix;
    vector<int> intervals;
    int size;
    int subsize;
    vector<bool> atomsInSubregion;

    // Basic Constructor:
    TConnectionWithIntervals(TSprinkling *Sprinkling) {
        atomsInSubregion.resize(0);
        causalmatrix.resize(0);
        linkmatrix.resize(0);
        intervals.resize(0);
        size = 0;
        subsize = 0;
        createCausalMatrixFromSprinkling(Sprinkling);
    }

    /** The link matrix requires a lot of effort to be computed.
      * Thus for applications where it is not required it can be skipped by
      * calling with needLinks=false. */
    TConnectionWithIntervals(TSprinkling *Sprinkling, bool needLinks, bool needSubregion) {
        atomsInSubregion.resize(0);
        causalmatrix.resize(0);
        linkmatrix.resize(0);
        intervals.resize(0);
        size = 0;
        subsize = 0;
        createCausalMatrixFromSprinkling(Sprinkling);
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
    void createCausalMatrixFromSprinkling(TSprinkling *Sprinkling) {
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

    /* INTERVAL METHODS */

    /*  This counts the number of inclusive order intervals of cardinality
     *  n + 1 for given n. Special cases:
     *  n = 1 gives the number of elements
     *  n = 2 gives the number of links  */
    int countIntervalsOfSize(int n) {
        if (n == 0) throw "Argument must be non-zero";
        if (n == 1) return this->size;
        if (this->size == 0) return 0;
        if (n > this->size) return 0;
        if (this->size == 1) {
            if (n == 1) return 1;
            else return 0;
        }
        int i, j, k, currentint;
        int out = 0;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (causal(i, j)) {
                    currentint = 0;
                    k = i + 1;
                    while (k < j && currentint < n) {
                        if (causal(i, k) && causal(k, j)) currentint++;
                        k++;
                    }
                    if (currentint == n - 2) out++;
                }
            }
        }
        return out;
    }

    /*  Spits out a vector whose entries are n-vectors of integers
     *  that represent the inclusive order intervals of cardinality n */
    vector<vector<int> > getIntervalsOfSize(int n) {
        if (n < 2) throw "Argument must be greater than 1";
        vector<vector<int> > out;
        vector<int> interval(n);
        int c_length;
        int i, j, k;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (causal(i, j)) {
                    c_length = 0;
                    interval.clear();
                    interval.push_back(i);
                    k = i + 1;
                    while (k < j && c_length < n) {
                        if (causal(i, k) && causal(k, j)) {
                            interval.push_back(k);
                            c_length++;
                        }
                        k++;
                    }
                    if (c_length == n - 2) {
                        interval.push_back(j);
                        out.push_back(interval);
                    }
                }
            }
        }
        return out;
    }

    /*  Given a subregion X find card(N_1(X|X)-N_1(X|M)) */
    int countD1X() {
        int count = 0;
        int i, j, k;
        int inside, outside;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (inSubregion(i) && inSubregion(j) && causal(i, j)) {
                     inside = 0;
                    outside = 0;
                    k = i + 1;
                    while (k < j && inside == 0) {
                        if (causal(i, k) && causal(k, j)) {
                            if (inSubregion(k))  inside++;
                            else                outside++;
                        }
                        k++;
                    }
                    if (inside == 0 && outside > 0) count++;
                }
            }
        }
        return count;
    }

    /*  Given a subregion X get N_1(X|X)-N_1(X|M) and print out all
     *  the intervals in a vector of vectors */
    vector<vector<int> > getD1X() {
        vector<vector<int> > vec;
        vector<int> currentint;
        int i, j, k;
        int inside, outside;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (inSubregion(i) && inSubregion(j) && causal(i, j)) {
                    currentint.clear();
                    currentint.push_back(i);
                     inside = 0;
                    outside = 0;
                    k = i + 1;
                    while (k < j && inside == 0) {
                        if (causal(i, k) && causal(k, j)) {
                            currentint.push_back(k);
                            if (inSubregion(k))  inside++;
                            else                outside++;
                        }
                        k++;
                    }
                    currentint.push_back(j);
                    if (inside == 0 && outside > 0) vec.push_back(currentint);
                }
            }
        }
        return vec;
    }

    /*  Given a subregion X and Y = M\X find card(D1(Y)) where
     *  D1(Y) = N_1(Y|Y)-N_1(Y|M) */
    int countD1Y() {
        int count = 0;
        int i, j, k;
        int inside, outside;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (!inSubregion(i) && !inSubregion(j) && causal(i, j)) {
                     inside = 0;
                    outside = 0;
                    k = i + 1;
                    while (k < j && inside == 0) {
                        if (causal(i, k) && causal(k, j)) {
                            if (!inSubregion(k))  inside++;
                            else                outside++;
                        }
                        k++;
                    }
                    if (inside == 0 && outside > 0) count++;
                }
            }
        }
        return count;
    }

    /*  Given a subregion X and Y=M\X get N_1(Y|Y)-N_1(Y|M) and print out all
     *  the intervals in a vector of vectors */
    vector<vector<int> > getN1Y() {
        vector<vector<int> > vec;
        vector<int> currentint;
        int i, j, k;
        int inside, outside;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (!inSubregion(i) && !inSubregion(j) && causal(i, j)) {
                    currentint.clear();
                    currentint.push_back(i);
                     inside = 0;
                    outside = 0;
                    k = i + 1;
                    while (k < j && inside == 0) {
                        if (causal(i, k) && causal(k, j)) {
                            currentint.push_back(k);
                            if (!inSubregion(k))  inside++;
                            else                outside++;
                        }
                        k++;
                    }
                    currentint.push_back(j);
                    if (inside == 0 && outside > 0) vec.push_back(currentint);
                }
            }
        }
        return vec;
    }

    int countIntervalsXXOfSize(int n) {
        if (n == 0) throw "Argument must be non-zero";
        if (n == 1) return this->subsize;
        if (this->subsize == 0) return 0;
        if (n > this->subsize) return 0;
        if (this->subsize == 1) {
            if (n == 1) return 1;
            else return 0;
        }
        int i, j, k, currentint;
        int out = 0;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (causal(i, j) && inSubregion(i) && inSubregion(j)) {
                    currentint = 0;
                    k = i + 1;
                    while (k < j && currentint < n) {
                        if (inSubregion(k) && causal(i, k) && causal(k, j)) currentint++;
                        k++;
                    }
                    if (currentint == n - 2) out++;
                }
            }
        }
        return out;
    }

    int countIntervalsXMOfSize(int n) {
        if (n == 0) throw "Argument must be non-zero";
        if (n == 1) return this->subsize;
        if (this->subsize == 0) return 0;
        if (n > this->subsize) return 0;
        if (this->subsize == 1) {
            if (n == 1) return 1;
            else return 0;
        }
        int i, j, k, currentint;
        int out = 0;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (causal(i, j) && inSubregion(i) && inSubregion(j)) {
                    currentint = 0;
                    k = i + 1;
                    while (k < j && currentint < n) {
                        if (causal(i, k) && causal(k, j)) currentint++;
                        k++;
                    }
                    if (currentint == n - 2) out++;
                }
            }
        }
        return out;
    }

    int countIntervalsYYOfSize(int n) {
        int ysize = this->size - this->subsize;
        if (n == 0) throw "Argument must be non-zero";
        if (n == 1) return this->size-this->subsize;
        if (ysize == 0) return 0;
        if (n > ysize) return 0;
        if (ysize == 1) {
            if (n == 1) return 1;
            else return 0;
        }
        int i, j, k, currentint;
        int out = 0;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (causal(i, j) && !inSubregion(i) && !inSubregion(j)) {
                    currentint = 0;
                    k = i + 1;
                    while (k < j && currentint < n) {
                        if (!inSubregion(k) && causal(i, k) && causal(k, j)) currentint++;
                        k++;
                    }
                    if (currentint == n - 2) out++;
                }
            }
        }
        return out;
    }

    int countIntervalsYMOfSize(int n) {
        int ysize = this->size - this->subsize;
        if (n == 0) throw "Argument must be non-zero";
        if (n == 1) return this->subsize;
        if (ysize == 0) return 0;
        if (n > ysize) return 0;
        if (ysize == 1) {
            if (n == 1) return 1;
            else return 0;
        }
        int i, j, k, currentint;
        int out = 0;
        for (i = 0; i < size - 1; i++) {
            for (j = i + 1; j < size; j++) {
                if (causal(i, j) && !inSubregion(i) && !inSubregion(j)) {
                    currentint = 0;
                    k = i + 1;
                    while (k < j && currentint < n) {
                        if (causal(i, k) && causal(k, j)) currentint++;
                        k++;
                    }
                    if (currentint == n - 2) out++;
                }
            }
        }
        return out;
    }

    /* computes the cardinality of the (inclusive) causal interval between
     * elements i and j (the endpoints i and j are counted). */
    int cardinality(int i, int j) {
        int result = 0;
        int a;
        for (a = i + 1; a < j; a++) {
            if (causal(i, a) && causal(a, j))
            result++;
        }
        return result + 2;
    }

    //TODO: finish this sheeeeeyat! and use dp like a baws.
//    int distance(int i, int j) {
//        if(!causal(i, j)) throw  "i must precede j";
//        out = 0;
//        for(int k = i + 1; k < j; k++) {
//
//        }
//        return out;
//    }

    // Counts the number of inclusive order intervals of cardinality n + 1 for all n.
    void storeAllIntervals() {
        this->intervals.assign(this->size, 0);
        this->intervals[0] = this->size;
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                if (causal(i, j)) this->intervals[cardinality(i, j) - 1]++;
            }
        }
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
    void writeCausal(string filename) {
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

    // Same for link matrix.
    void writeLink(ofstream *os) {
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
