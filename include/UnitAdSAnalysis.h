#pragma once

#include <fstream>
#include <memory>
#include <new>
#include <string>
#include <vector>

#include "Tools.h"
#include "UnitConnection.h"
#include "UnitGeometry_AntiDeSitter2.h"
#include "UnitGeometry.h"
#include "UnitSprinkling.h"

using namespace std;

class TAdSAnalysis {

private:
    TSprinkling* Sprinkling;
    vector<TPoint*> vertex;

public:
    TAdSAnalysis(TSprinkling* _Sprinkling) {
//        if (Sprinkling->geometry != TGeometry_AntiDeSitter2) {
//            throw "Can only use TAdSAnalysis class with AdS2 geometry!";
//        }
        this->Sprinkling = _Sprinkling;
        this->vertex = Sprinkling->points;
    }

    // empty constructor for analysis();
    TAdSAnalysis() {
        Sprinkling = NULL;
    };

    // does a precede b after projection onto the boundary?
    // since the projection just puts the x-value of an event to
    // plus or minus pi/2 depending on the sign of x, a point a
    // will precede a point b on the boundary iff t(b) > t(a)
    // and their x-values have the same sign in the bulk.
    bool precOnBoundary(TPoint *_a, TPoint *_b) {
        Tp_AntiDeSitter2 *a = (Tp_AntiDeSitter2*) _a;
        Tp_AntiDeSitter2 *b = (Tp_AntiDeSitter2*) _b;
        if ((b->t > a->t) && ((b->x >= 0) ^ (a->x < 0))) return true;
        else return false;
    }

    // does a precede b after projection onto the boundary? see above.
    bool precOnBoundary(int i, int j) {
        Tp_AntiDeSitter2 *a = (Tp_AntiDeSitter2*) Sprinkling->points[i];
        Tp_AntiDeSitter2 *b = (Tp_AntiDeSitter2*) Sprinkling->points[j];
        if ((b->t > a->t) && ((b->x >= 0) ^ (a->x < 0))) return true;
        else return false;
    }

    // print causal matrix after projection onto the boundary to file
    void writeCausalOnBoundary(string filename) {
        ofstream os;
        os.open(filename.c_str());
        for (int i = 0; i < Sprinkling->size; i++) {
            for (int j = 0; j < Sprinkling->size; j++) {
                os << " " << this->precOnBoundary(i, j);
            }
            os << endl;
        }
        cout << "Printed causal matrix to the file \"" << filename << "\"."
                << endl;
        os.close();
    }

    // evaluate fraction of relations that are preserverd after
    // projection onto the boundary
    double fractionOfPreservedRelations() {
        double preservedrels = 0.0;
        double totalrels = 0.0;
        for (int i = 0; i < Sprinkling->size; i++) {
            for (int j = 0; j < Sprinkling->size; j++) {
                if (Sprinkling->prec(i, j)) {
                    totalrels++;
                    if (this->precOnBoundary(i, j)) preservedrels++;
                }
            }
        }
        return preservedrels / totalrels;
    }

    void analysis() {
        double rmax;
        TGeometry *Geometry;
        TSprinkling *Sprinkling;
        TConnection *Connection;
        TAdSAnalysis *AdSAnalysis;

        double epsilon, rhobulk, rhobound, Vbound;
        ofstream os;
        os.open("analysis.dat");
        os << "rmax \t size \t rhobulk \t rhobound \t PR" << endl;
        for (int it = 0; it < 100; it++) {
            epsilon = 0.5 - 0.5 * (double) it / 100;
            rmax = PI / 2.0 - epsilon - 0.1;
            Vbound = 2.0/cos(rmax);

            Geometry = new TGeometry_AntiDeSitter2(rmax, 0.0,1.0);
            Sprinkling = new TSprinkling(1000.0*Vbound, false, Geometry);
            Connection = new TConnection(Sprinkling);
            AdSAnalysis = new TAdSAnalysis(Sprinkling);

            //AdSAnalysis->writeCausalOnBoundary("141013ads2-boundary-cmatrix.dat");

            cout << "step " << it << "\t rmax = " << rmax << "\t N = " << Sprinkling->size << endl;
            rhobulk = Sprinkling->size/Geometry->V;
            rhobound = Sprinkling->size/Vbound;
            os << rmax << "\t";
            os << Sprinkling->size << "\t";
            os << rhobulk << "\t";
            os << rhobound << "\t";
            os << AdSAnalysis->fractionOfPreservedRelations() << endl;
        };
    }

};


