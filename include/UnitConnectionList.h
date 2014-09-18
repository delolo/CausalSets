/*----------------------------------------------------------------------------
 *
 *  Written:       25/05/2014
 *  Last updated:  25/05/2014
 *
 *  Same as TConnection but with an adjacency list implementation instead
 *  of an adjacency matrix implementation.
 *
 *---------------------------------------------------------------------------*/

#pragma once

//#include <__tree>
//#include <cmath>
#include <fstream>
//#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

//#include "Tools.h"
#include "UnitSprinkling.h"

using namespace std;

class TCauset {
private:
    vector<vector<int> > futurelist;
    vector<vector<int> > futurelinklist;
    int size;

public:
    // Basic Constructor:
    TCauset(TSprinkling *Sprinkling) {
        size = 0;
    }

    void causalFromSprinkling(TSprinkling *Sprinkling) {
        this->size = Sprinkling->size;
        int i, j;
        futurelist.resize(size);
        for (i = 0; i < size; i++) {
            futurelist[i].resize(0);
            for (j = 0; j < size - i - 1; j++) {
                if(Sprinkling->prec(i, j)) futurelist[i].push_back(j);
            }
        }
    }

    /* BASIC CAUSAL RELATIONS METHODS */

public:
    /*  For applications it is much more convenient to have something that looks more
     *  like a full matrix. For this the function causal(i,j) can be used. causal(i,j)
     *  returns whether element i in the Sprinkling->vertex is to the causal past
     *  of element j. */
    bool causal(int i, int j) {
        if (j <= i) return false;
        return find(futurelist[i].begin(), futurelist[i].end(), j)!=futurelist[i].end();
    }

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

    // To decide whether two causally related elements are linked one must go over all
    // potential elements in between and rule them out. Cost is O(N^3).
    void linkFromSprinkling(TSprinkling* Sprinkling) {}

    /*  See causal(i,j). This is the same for links. */
    bool link(int i, int j) {
        return false;
    }

    /* IO STUFF */

    // adjacency list as a string
    string toString()
    {
        std::ostringstream str;
        str << "Causal Matrix:\n";
        for (int i = 0; i < size; i++) {
            str << "\n";
            for (vector<int>::iterator iter = futurelist[i].begin(); iter != futurelist[i].end(); iter++) {
                str << *iter << " ";
            }
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
}
;

// THomotopyMatrix is literally identical to TConnection without the link part.
// Except that the usual causal relation is replaced by the relation that counts
// the number of equivalence classes of causal paths. This can be used to compute
// for example the propagator on a cylinder.
