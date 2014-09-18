/*----------------------------------------------------------------
 *
 * 	Written:       10/10/2013
 *  Last updated:  10/10/2013
 *
 *  Sample main file.
 *
 *----------------------------------------------------------------*/

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "../include/Statistics.h"
#include "../include/Timer.h"
#include "../include/Tools.h"
#include "../include/UnitAdSAnalysis.h"
#include "../include/UnitConnection.h"
#include "../include/UnitConnectionList.h"
#include "../include/UnitConnectionWithIntervals.h"
#include "../include/UnitExtrinsicCurvature.h"
#include "../include/UnitExtrinsicCurvature2.h"
#include "../include/UnitExtrinsicCurvatureSmeared.h"
#include "../include/UnitGeometry_AntiDeSitter2.h"
#include "../include/UnitGeometry_AntiDeSitter3.h"
#include "../include/UnitGeometry_dS2ConformalSlab.h"
#include "../include/UnitGeometry_Minkowski2_Rectangle.h"
#include "../include/UnitGeometry_MinkowskiD_Rectangle.h"
#include "../include/UnitGeometry_Minkowski2_Diamond.h"
#include "../include/UnitGeometry_Minkowski3_Cylinder.h"
#include "../include/UnitGeometry_EinsteinStaticUniverseD.h"
#include "../include/UnitGeometry_CFTimeSquaredD_Rectangle.h"
#include "../include/UnitGeometry_SpatialLinear2_Rectangle.h"
#include "../include/UnitGeometry_SpatialExponential2_Rectangle.h"
#include "../include/UnitGeometry_dS4Global.h"
#include "../include/UnitGeometry.h"
#include "../include/UnitPropagator.h"
#include "../include/UnitSprinklingLabeled.h"
#include "../include/UnitSprinkling.h"
#include "../include/UnitPropagators.h"

using namespace std;

void testPropagator() {

    using namespace std;

    cout << "alpha=1.0\n";
    double t0 = 0.0;
    double t1 = PI / 2.0 - 0.1;
    printf("t0 = %g, t1 = %g\n", t0, t1);
    int N;

    cout << "<N> of the Causal Set:        ";
    cin >> N;

    int s;
    TPoint *atom;
    TGeometry *Geometry;
    TSprinkling *Sprinkling;
    TConnectionWithIntervals *Connection;
    TPropagator *Propagator;
    Geometry = new TGeometry_Minkowski2_Rectangle;
    Sprinkling = new TSprinkling(N, Geometry);
    Connection = new TConnectionWithIntervals(Sprinkling);
    Propagator = new TPropagator(Sprinkling, 0.3);
    s = Sprinkling->size;

//write causal set coordinates and additional information into file
    cout << (Sprinkling->size) << endl;

    ofstream os;
    os.open("coords.txt");
    for (int n = 0; n < Sprinkling->size; n++) {
        atom = Sprinkling->points[n];
        os << atom->t << "\t" << ((Tp_Minkowski2_Rectangle*) atom)->x
                << endl;
    }
    os.close();
    cout << "printed coordinates" << endl;

//write information for the propagator algorithm into file

    ofstream data;
    data.open("data.txt");
    data << Geometry->V << endl;
    data.close();

//write causal matrix into file:

    Connection->writeCausal("causal.txt");

    double ** kret = Propagator->getRetardedPropagator();

    for (int i = 0; i < s; i++) {
        for (int j = 0; j < s; j++) {
            cout << std::abs(kret[i][j]) << " ";
        }
        cout << endl;
    }

//display
//
//	cout << "Number of causal set elements:\t" << Sprinkling->size << endl;
//	cout << "Volume of the causal diamond: \t" << Geometry->V << endl;
//	cout << "Density of elements: \t" << Sprinkling->size / Geometry->V << endl;

//	Geodesic Distances:

//	Clear memory to avoid segmentation faults.

    delete Geometry;
    delete Sprinkling;
//	delete Connection;
    delete Propagator;
}

void testSpherical() {
    int N, d;
    cout << "dimension d of the d-sphere?" << endl;
    cin >> d;
    cout << "how many points?" << endl;
    cin >> N;
    ofstream f;
    f.open("newsphere.txt");
    vector<double> coords;
    for (int i = 0; i < N; i++) {
        coords = rndSpherical(d, 1.0);
        cout << vecToString(coords) << endl;
        for (int k = 0; k < d + 1; k++) {
            f << coords[k] << " ";
        }
        f << "\n";
    }
    f.close();
}

void testHemispherical() {
    int N, d;
    cout << "dimension d of the d-sphere?" << endl;
    cin >> d;
    cout << "how many points?" << endl;
    cin >> N;
    ofstream f;
    f.open("newsphere.txt");
    vector<double> coords;
    for (int i = 0; i < N; i++) {
        coords = rndSpherical(d, 3.0);
        if (coords[0] < 0) coords[0] = -coords[0]; // project lower onto upper hemisphere
        for (int k = 0; k < d + 1; k++) {
            f << coords[k] << " ";
        }
        f << "\n";
    }
    f.close();
}

void testAntiDeSitter2() {
    int N;
    cout << "number of elements?";
    cin >> N;
    TGeometry *Geometry;
    TSprinkling *Sprinkling;

    double epsilon;
    ofstream os;
    os.open("analysis.dat");
    for (int it = 0; it <= 100; it++) {
        epsilon = 0.5 - 0.5 * (double) it / 100;
        cout << "epsilon = " << epsilon << endl;

        Geometry = new TGeometry_AntiDeSitter2(PI / 2.0 - epsilon, 0.0,
                1.0);

        Sprinkling = new TSprinkling(N, Geometry);
        TConnectionWithIntervals* Connection = new TConnectionWithIntervals(
                Sprinkling);

//    Sprinkling->writeToFile("141013ads2-coords.dat");
//    Connection->writeCausal("141013ads2-cmatrix.dat");

        TAdSAnalysis* AdSAnalysis = new TAdSAnalysis(Sprinkling);

//        AdSAnalysis->writeCausalOnBoundary("141013ads2-boundary-cmatrix.dat");

        os << PI / 2.0 - epsilon << "\t"
                << AdSAnalysis->fractionOfPreservedRelations() << endl;

        delete Connection;
        delete Geometry;
        delete Sprinkling;
    }
    os.close();

}

void testAntiDeSitter3() {
    int N;
    cout << "number of elements?" << endl;
    cin >> N;
    TGeometry* Geometry;
    TSprinkling* Sprinkling;
    Geometry = new TGeometry_AntiDeSitter3(PI / 2.0 - 0.15, 0.0, 1.0);

    Timer* timer1 = new Timer("sprinkling+connection");
    Sprinkling = new TSprinkling(N, Geometry);
    TConnectionWithIntervals* Connection = new TConnectionWithIntervals(
            Sprinkling);
    timer1->stop();

    Sprinkling->writeToFile("111013ads3-coords.dat");
    Connection->writeCausal("111013ads3-cmatrix.dat");

    delete Connection;
    delete Geometry;
    delete Sprinkling;

}

void testAdSAnalysis() {
    TAdSAnalysis* AdSA = new TAdSAnalysis();
    AdSA->analysis();
}

void testExtrinsicCurvatureDS2SingleRho() {
    // let's choose some time ts that defines the surface and choose two slabs to the
    // future and past of t=ts. let the slab above go from ts to ttop and the slab
    // below go from ts to tbottom. let's choose the two slabs such that they have equal
    // volumes, which requires tbottom = 1/(2/ts - ttop);

    // in this case K for any constant t surface is just 1/alpha. and the analytic
    // calculation gives minminusmax = -4/3 width/ts.

    double ttop = 2.0;
    double ts = 1.0;
    double tbottom = 1.0 / 1.5;
    double width = 2.0;
    cout << "starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_dS2ConformalSlab(1.0, ts, ttop,
            width);
    TGeometry* GeometryBelow = new TGeometry_dS2ConformalSlab(1.0, tbottom,
            ts,
            width);
    double density = 500.0;
    cout << "expected size of causet   = " << density * GeometryAbove->V * 2
            << endl;
    cout << "number of causets sampled = " << 10000 << endl;
    Timer* timer = new Timer("stats");
    TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(2, GeometryAbove,
            GeometryBelow, density);
    vector<double> stats = Curvature->getNewStats(10000);
    printf("MinMinusMax      = %4.4f +- %4.4f\n", stats[0], stats[2]);
    cout << "Analytic result  = " << -4.0 / 3.0 * width / ts << endl;
    timer->stop();
}

void testExtrinsicCurvatureDS2() {
    // let's choose some time ts that defines the surface and choose two slabs to the
    // future and past of t=ts. let the slab above go from ts to ttop and the slab
    // below go from ts to tbottom. let's choose the two slabs such that they have equal
    // volumes, which requires tbottom = 1/(2/ts - ttop);

    // in this case K for any constant t surface is just 1/alpha. and the analytic
    // calculation gives minminusmax = -4/3 width/ts.

    double ttop = 2.0;
    double ts = 1.0;
    double tbottom = 1.0 / 1.5;
    double width = 2.0;
    cout << "starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_dS2ConformalSlab(1.0, ts, ttop,
            width);
    TGeometry* GeometryBelow = new TGeometry_dS2ConformalSlab(1.0, tbottom,
            ts,
            width);
    double density;
    ofstream datafile;
    ofstream statfile;
    datafile.open("minminusmax-desitter-data.txt", ios::app);
    statfile.open("minminusmax-desitter-stat.txt", ios::app);
    for (int k = 14; k < 15; k++) {
        density = pow(2.0, (double) (k + 3));
        cout << "expected size of causet   = "
                << density * GeometryAbove->V * 2
                << endl;
        cout << "number of causets sampled = " << 1000 << endl;
        Timer* timer = new Timer("stats");
        TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(2,
                GeometryAbove,
                GeometryBelow, density);
        Curvature->printNewDataAndStats(1000, &datafile, &statfile);
        cout << "Analytic result  = " << -4.0 / 3.0 * width / ts << endl;
        timer->stop();
    }
    datafile.close();
    statfile.close();
}

void testExtrinsicCurvatureMink2() {
    cout << "Starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_Minkowski2_Rectangle(0.5, 1.0);
    TGeometry* GeometryBelow = new TGeometry_Minkowski2_Rectangle(0.5, 1.0);
    double density;
    ofstream datafile;
    ofstream statfile;
    datafile.open("140528minminusmax-mink2-data.txt", ios::app);
    statfile.open("140528minminusmax-mink2-stat.txt", ios::app);
    int k;
    for (k = 2; k < 35; k++) {
        density = pow(2.0, (double) k / 2.0);
        cout << "expected size of causet   = " << density * GeometryAbove->V * 2
                << endl;
        cout << "number of causets sampled = " << 100 << endl;
        Timer* timer = new Timer("stats");
        TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(2,
                GeometryAbove, GeometryBelow, density);
        Curvature->printNewDataAndStats(100, &datafile, &statfile);
        cout << "Analytic result  = " << 0.0 << endl;
        timer->stop();
        delete Curvature;
        delete timer;
    }
    datafile.close();
    statfile.close();
    delete GeometryAbove;
    delete GeometryBelow;
}

void testExtrinsicCurvatureMinkowskiDRectangle() {
    cout << "Starting..." << endl;

    int dimension = 4;
    int RUNS = 100;

    double L = 1.0;
    vector<double> xmax;
    xmax.resize(dimension - 1, L);
    cout << vecToString(xmax);

    TGeometry* GeometryAbove =
            new TGeometry_MinkowskiD_Rectangle(dimension, xmax, 0.5);
    TGeometry* GeometryBelow =
            new TGeometry_MinkowskiD_Rectangle(dimension, xmax, 0.5);
    ofstream datafile, statfile;
    stringstream datafilename, statfilename;
    datafilename << currentDate()
            << "minminusmax-mink"
            << dimension
            << "-runs"
            << RUNS
            << "-L"
            << fixed << setprecision(0) << L
            << "-T1-data.txt";
    statfilename << currentDate()
            << "minminusmax-mink"
            << dimension
            << "-runs"
            << RUNS
            << "-L"
            << fixed << setprecision(0) << L
            << "-T1-stat.txt";
    datafile.open(datafilename.str(), ios::app);
    statfile.open(statfilename.str(), ios::app);

    double nexp;
    double rhoexp;
    for (int k = 10; k < 11; k++) {
        nexp = pow(2.0, (double) k / 2.0);
        rhoexp = nexp / (GeometryAbove->V + GeometryBelow->V);
        cout << "expected density          = " << rhoexp << endl;
        cout << "number of causets sampled = " << RUNS << endl;

        Timer* timer = new Timer("stats");
        TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(dimension,
                GeometryAbove, GeometryBelow, rhoexp);
        Curvature->printNewDataAndStats(RUNS, &datafile, &statfile);
        cout << "Analytic result  = " << 0.0 << endl;
        timer->stop();

        delete Curvature;
        delete timer;
    }
    datafile.close();
    statfile.close();
    delete GeometryAbove;
    delete GeometryBelow;
}

void testExtrinsicCurvatureMinkowski3Rectangle_DependenceOnL() {
    cout << "Starting..." << endl;

    int dimension = 3;

    stringstream datafilename, statfilename;
    datafilename << currentDateTime()
            << "minminusmax-mink"
            << dimension
            << "N16K-Lvar-data.txt";
    statfilename << currentDateTime()
            << "minminusmax-mink"
            << dimension
            << "N16K-Lvar-stat.txt";

    for (int l = 0; l < 6; l++) {
        vector<double> xmax;
        xmax.resize(dimension - 1, 1 + 10.0 * l);
        cout << vecToString(xmax);

        TGeometry* GeometryAbove =
                new TGeometry_MinkowskiD_Rectangle(dimension, xmax, 0.5);
        TGeometry* GeometryBelow =
                new TGeometry_MinkowskiD_Rectangle(dimension, xmax, 0.5);
        ofstream datafile;
        ofstream statfile;
        datafile.open(datafilename.str(), ios::app);
        statfile.open(statfilename.str(), ios::app);

        double nexp;
        double rhoexp;
        int k;
        for (k = 28; k < 29; k++) {
            nexp = pow(2.0, (double) k / 2.0);
            rhoexp = nexp / (GeometryAbove->V + GeometryBelow->V);

            cout << "expected density          = " << rhoexp << endl;
            cout << "number of causets sampled = " << 100 << endl;

            Timer* timer = new Timer("stats");
            TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(dimension,
                    GeometryAbove, GeometryBelow, rhoexp);
            Curvature->printNewDataAndStats(100, &datafile, &statfile);
            cout << "Analytic result  = " << 0.0 << endl;
            timer->stop();

            delete Curvature;
            delete timer;
        }
        datafile.close();
        statfile.close();
        delete GeometryAbove;
        delete GeometryBelow;
    }
}

void testExtrinsicCurvatureTimeSquaredDRectangle() {
    cout << "Starting..." << endl;

    int dimension = 4;

    double L = 1.0;
    vector<double> xmax;
    xmax.resize(dimension - 1, L);

    ofstream datafile, statfile;
    stringstream datafilename, statfilename;
    datafilename << currentDate()
            << "minminusmax-tlinear"
            << dimension
            << "-L"
            << fixed << setprecision(0) << L
            << "-P"
            << "12"
            << "-ts1TO20-data.txt";
    statfilename << currentDate()
            << "minminusmax-tlinear"
            << dimension
            << "-L"
            << fixed << setprecision(0) << L
            << "-P"
            << "14"
            << "-ts1TO20-stat.txt";
    datafile.open(datafilename.str(), ios::app);
    statfile.open(statfilename.str(), ios::app);

    datafile << "# This is data for an N = 2^12 sprinkling into the"  << endl
             << "# [0,1]^3 x [1,20] slab of the TimeLinear4 spacetime"<< endl
             << "# for ts=4,...,17." << endl;
    statfile << "# This is data for an N = 2^12 sprinkling into the"  << endl
             << "# [0,1]^3 x [1,20] slab of the TimeLinear4 spacetime"<< endl
             << "# for ts=4,...,17." << endl;

    for (int t = 4; t < 17; t++) {
        double T = (double) t;
        TGeometry* GeometryBelow = new TGeometry_CFTimeSquaredD_Rectangle(
                dimension, xmax, 1.0, T);
        TGeometry* GeometryAbove = new TGeometry_CFTimeSquaredD_Rectangle(
                dimension, xmax, T, 20.0);
        double nexp;
        double rhoexp;
        for (int k = 28; k < 29; k++) {
            nexp = pow(2.0, (double) k / 2.0);
            rhoexp = nexp / (GeometryAbove->V + GeometryBelow->V);

            cout << "expected density          = " << rhoexp << endl;
            cout << "number of causets sampled = " << 100 << endl;

            Timer* timer = new Timer("stats");
            TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(dimension,
                    GeometryAbove, GeometryBelow, rhoexp);
            statfile << t << "\t";
            datafile << t << "\t";
            Curvature->printNewDataAndStats(100, &datafile, &statfile);
            cout << "Analytic result  = " << 0.0 << endl;
            timer->stop();

            delete Curvature;
            delete timer;
        }
        delete GeometryAbove;
        delete GeometryBelow;
    }
    datafile.close();
    statfile.close();
}

void debugExtrinsicCurvatureSmearedMink2() {
    cout << "Starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_Minkowski2_Rectangle(0.5, 1.0);
    TGeometry* GeometryBelow = new TGeometry_Minkowski2_Rectangle(0.5, 1.0);
    double density;
    double epsilon;
//    ofstream datafile;
//    ofstream statfile;
//    datafile.open("140526minminusmax-test-data.txt", ios::app);
//    statfile.open("140526minminusmax-test-stat.txt", ios::app);
    density = pow(2.0, 8.0);
    epsilon = 1.0 / 25.0;
    cout << "expected size of causet   = " << density * GeometryAbove->V * 2
            << endl;
    cout << "number of causets sampled = " << 1 << endl;
    Timer* timer = new Timer("stats");
    TExtrinsicCurvatureSmeared* Curvature =
            new TExtrinsicCurvatureSmeared(2, GeometryAbove, GeometryBelow,
                    density);
    cout << Curvature->getNewMinMinusMax(epsilon) << endl;
    cout << "Analytic result  = " << 0.0 << endl;
    timer->stop();
    delete Curvature;
    delete timer;
}

void testExtrinsicCurvatureSmearedMink2() {
    cout << "Starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_Minkowski2_Rectangle(0.5, 1.0);
    TGeometry* GeometryBelow = new TGeometry_Minkowski2_Rectangle(0.5, 1.0);
    double density;
    double epsilon;
    int RUNS = 100;
    ofstream datafile;
    ofstream statfile;
    datafile.open("140528minminusmax-smeared-lk1-data.txt", ios::app);
    statfile.open("140528minminusmax-smeared-lk1-stat.txt", ios::app);
    int k;
    for (k = 8; k < 34; k++) {
        density = pow(2.0, (double) k / 2.0);
        epsilon = (1.0 / density) / 1.0;
        cout << "epsilon = " << epsilon << endl;
        cout << "expected size of causet   = "
                << density * GeometryAbove->V * 2 << endl;
        cout << "number of causets sampled = "
                << RUNS << endl;
        Timer* timer = new Timer("stats");
        TExtrinsicCurvatureSmeared* Curvature =
                new TExtrinsicCurvatureSmeared(2, GeometryAbove, GeometryBelow,
                        density);
        Curvature->printNewDataAndStats(RUNS, epsilon, &datafile, &statfile);
        cout << "Analytic result  = " << 0.0 << endl;
        timer->stop();
        delete Curvature;
        delete timer;
    }
    datafile.close();
    statfile.close();
    delete GeometryAbove;
    delete GeometryBelow;
}

void testExtrinsicCurvatureMink3Cylinder() {
    cout << "Starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_Minkowski3_Cylinder();
    TGeometry* GeometryBelow = new TGeometry_Minkowski3_Cylinder();
    double density;
    ofstream datafile;
    ofstream statfile;
    datafile.open("140529minminusmax-mink3-data.txt", ios::app);
    statfile.open("140529minminusmax-mink3-stat.txt", ios::app);
    int k;
    for (k = 4; k < 35; k++) {
        density = pow(2.0, (double) k / 2.0);
        cout << "expected size of causet   = " << density * GeometryAbove->V * 2
                << endl;
        cout << "number of causets sampled = " << 100 << endl;
        Timer* timer = new Timer("stats");
        TExtrinsicCurvature* Curvature = new TExtrinsicCurvature(3,
                GeometryAbove, GeometryBelow, density);
        Curvature->printNewDataAndStats(100, &datafile, &statfile);
        cout << "Analytic result  = " << 0.0 << endl;
        timer->stop();
        delete Curvature;
        delete timer;
    }
    datafile.close();
    statfile.close();
    delete GeometryAbove;
    delete GeometryBelow;
}

void testExtrinsicCurvatureSmearedMink3() {
    cout << "Starting..." << endl;
    TGeometry* GeometryAbove = new TGeometry_Minkowski3_Cylinder();
    TGeometry* GeometryBelow = new TGeometry_Minkowski3_Cylinder();
    double density;
    double epsilon;
    int RUNS = 100;
    ofstream datafile;
    ofstream statfile;
    datafile.open("140529minminusmax-cyl3-smeared-lk1-data.txt", ios::app);
    statfile.open("140529minminusmax-cyl3-smeared-lk1-stat.txt", ios::app);
    int k;
    for (k = 8; k < 28; k++) {
        density = pow(2.0, (double) k / 2.0);
        epsilon = 1.0 / (density * 9.0);
        cout << "epsilon = " << epsilon << endl;
        cout << "expected size of causet   = "
                << density * GeometryAbove->V * 2 << endl;
        cout << "number of causets sampled = "
                << RUNS << endl;
        Timer* timer = new Timer("stats");
        TExtrinsicCurvatureSmeared* Curvature =
                new TExtrinsicCurvatureSmeared(3, GeometryAbove, GeometryBelow,
                        density);
        Curvature->printNewDataAndStats(RUNS, epsilon, &datafile, &statfile);
        cout << "Analytic result  = " << 0.0 << endl;
        timer->stop();
        delete Curvature;
        delete timer;
    }
    datafile.close();
    statfile.close();
    delete GeometryAbove;
    delete GeometryBelow;
}

void testExtrinsicCurvature2MinkowskiDRectangle() {
    cout << "Starting..." << endl;

    int dimension = 4;
    int RUNS = 100;

    double L = 1.0;
    vector<double> xmax;
    xmax.resize(dimension - 1, L);

    TGeometry* Geometry = new TGeometry_MinkowskiD_Rectangle(dimension, xmax, 1.0);
    ofstream datafile, statfile;
    stringstream datafilename, statfilename;
    datafilename << currentDate()
            << "schnitzel-mink"
            << dimension
            << "-runs"
            << RUNS
            << "-L"
            << fixed << setprecision(0) << L
            << "-T1-data.txt";
    statfilename << currentDate()
            << "schnitzel-mink"
            << dimension
            << "-runs"
            << RUNS
            << "-L"
            << fixed << setprecision(0) << L
            << "-T1-stat.txt";
    datafile.open(datafilename.str().c_str(), ios::app);
    statfile.open(statfilename.str().c_str(), ios::app);

    double nexp;
    double rhoexp;
    for (int k = 10; k < 11; k++) {
        nexp = pow(2.0, (double) k / 2.0);
        rhoexp = nexp / Geometry->V;
        cout << "expected density          = " << rhoexp << endl;
        cout << "number of causets sampled = " << RUNS << endl;

        Timer* timer = new Timer("stats");
        TExtrinsicCurvature2* Curvature =
                new TExtrinsicCurvature2(dimension, Geometry, rhoexp);
        Curvature->printNewDataAndStats(RUNS, &datafile, &statfile);
        cout << "Analytic result  = " << 0.0 << endl;
        timer->stop();

        delete Curvature;
        delete timer;
    }
    datafile.close();
    statfile.close();
    delete Geometry;
}

void debugSprinkling() {
    double ttop = 2.0;
    double ts = 1.0;
    double width = 2.0;
    int rho = 4;
    TGeometry* GeometryAbove = new TGeometry_dS2ConformalSlab(1.0, ts, ttop,
            width);
    for (int k = 0; k < 10; k++) {
        TSprinkling* Sprinkling = new TSprinkling(
                rho * GeometryAbove->V,
                GeometryAbove);
        //cout << "successfully sprinkled " << Sprinkling->size << " elements" << endl;
        delete Sprinkling;
    }
}

/* Generic test for sprinklings */
void testSprinkling() {

    int dimension = 3;

    vector<double> xmax;
    xmax.resize(dimension - 1, 5.0);

    TGeometry* Geometry = new TGeometry_MinkowskiD_Rectangle(3);
    cout << "volume = " << Geometry->V << endl;
    TSprinklingLabeled* Sprinkling;

    for (int RUN = 0; RUN < 1; RUN++) {
        Sprinkling = new TSprinklingLabeled(10, false, Geometry);
        Sprinkling->writeToFile("newtest.txt");

        cout << Sprinkling->toString();
        cout << endl << "ERASING A POINT" << endl;
        Sprinkling->remove(2);
        Sprinkling->remove(3);
        cout << Sprinkling->toString() << endl << endl;
        Sprinkling->reLabel();
        cout << Sprinkling->toString() << endl << endl;
        cout << Sprinkling->toString();

        delete Sprinkling;
    }
}

/* Generic test for sprinklings and causal matrices */
void testConnection() {
    int dimension = 2;
    TGeometry* Geometry = new TGeometry_EinsteinStaticUniverseD(dimension, 1.0, 0.0, 10.0);
    TSprinkling* Sprinkling;
    TConnection* Connection;

    for (int RUN = 0; RUN < 1; RUN++) {
        Sprinkling = new TSprinkling(50, false, Geometry);
        Sprinkling->writeToFile("testsprinkling.dat");
        Connection = new TConnection(Sprinkling, false, false);
        Connection->writeCausalMatrix("testcausal.dat");
        cout << Connection->toString();
        delete Sprinkling;
        //delete Connection;
    }
}

void testSprinklingLabeled() {

    int dimension = 2;

    vector<double> xmax;
    xmax.resize(dimension - 1, 5.0);

    TGeometry* Geometry = new TGeometry_Minkowski2_Rectangle();
    cout << "volume = " << Geometry->V << endl;
    TSprinkling* Sprinkling;

    for (int RUN = 0; RUN < 1; RUN++) {
        Sprinkling = new TSprinkling(10, false, Geometry);
        Sprinkling->writeToFile("newtest.txt");
        delete Sprinkling;
    }
}

void testSpatialExponentialSprinkling() {

    double rho = 1.0;
    TGeometry* Geometry;
    TSprinkling* Sprinkling;

    Timer* timer = new Timer("sprinkling");

    for (int RUN = 1; RUN < 4; RUN++) {
        Geometry = new TGeometry_SpatialLinear2_Rectangle(
                1.0 / (double) RUN);
        //rho = pow(2.0, RUN);
        timer->reset();
//        cout << "i was at 2\n";
        Sprinkling = new TSprinkling(rho, Geometry);
        cout << "rho = " << rho << endl;
        cout << "N   = " << Sprinkling->size << endl;
        cout << "V   = " << Geometry->V << endl;
//        cout << "hi" << endl;
        timer->stop();
        timer->reset();
        //Connection = new TConnectionWithIntervals(Sprinkling);
        //Connection->writeCausal("thiscausal.txt");
        timer->stop();
    }
//    delete Connection;
    delete Sprinkling;
}

void testIntervalsInCylinder() {
    int N = pow(2, 8);
    double subradius = 0.5;
    TGeometry* Geometry = new TGeometry_Minkowski3_Cylinder(subradius);
    ofstream os, os2;
    os.open("intervalTest-data.dat");

    TSprinkling* Sprinkling = new TSprinkling(N, false, Geometry);
    Sprinkling->writeToFile("intervalTest-sprinkling.txt");
    TConnectionWithIntervals* Connection = new TConnectionWithIntervals(
            Sprinkling, false, true);

    /* Calculate fraction of points in cylinder*/

    cout << "FRACTION = " << (double) Connection->subsize / Connection->size
            << endl;
    cout << "N1 = " << Connection->countIntervalsOfSize(2) << endl;
    Timer* timer2 = new Timer("countOrderIntervals");
    //vector<vector<int> > ints = Connection->getIntervalsOfCardinality(20);
    //for(int i = 0; i < Connection->size; i++) {if(Connection->inSubregion(i)) cout << i << endl;};
    cout << "N1(INSIDE) = " << Connection->countD1X() << endl;
    timer2->stop();

    delete Connection;
    delete Sprinkling;
    delete Geometry;
}

void testD1_Minkowski3_Cylinder() {
    for (int power = 1; power < 14; power++) {
        int N = pow(2, power);
        cout << "N = " << N << endl;
        double subradius = 0.3;
        double rho;
        TGeometry* Geometry;
        TSprinkling* Sprinkling;
        TConnectionWithIntervals* Connection;
        Statistics* myStats = new Statistics("D1");

        Geometry = new TGeometry_Minkowski3_Cylinder(subradius);

        ofstream datastream("cylinder-D1-data-subradius=0.3.txt", ios::app);

        for (int RUN = 0; RUN < 10; RUN++) {
            Sprinkling = new TSprinkling(N, false, Geometry);
            Connection = new TConnectionWithIntervals(Sprinkling, false,
                    true);

            rho = Sprinkling->size / Geometry->V;
            double myD1 = (Connection->countD1X() - Connection->countD1Y())
                    * pow(rho, -1.0 / 3.0);
            myStats->push(myD1);
            cout << myStats->getTime();
            delete Connection;
            delete Sprinkling;
        }
        delete Geometry;
        cout << "STATS = " << vecToString(myStats->getStats()) << endl;
        datastream << 2 << "^" << power << "\t";
        datastream << rho << "\t";
        datastream << myStats->toString();
    }
}

void testD1_Minkowski2_Rectangle() {
    for (int power = 11; power < 17; power++) {
        int N = pow(2, power);
        cout << "N = " << N << endl;
        std::vector<double> partition(2);
        partition[0] = 0.5;
        TGeometry* Geometry;
        TSprinkling* Sprinkling;
        TConnectionWithIntervals* Connection;
        Statistics* myStats = new Statistics("D1");
        Statistics* rhoStats = new Statistics("rhoStats");

        cout << "HERE" << endl;

        Geometry = new TGeometry_Minkowski2_Rectangle(1, partition);

        ofstream datastream("rectangle-D1-newdata.txt", ios::app);

        for (int RUN = 0; RUN < 50; RUN++) {
            Sprinkling = new TSprinkling(N, Geometry);
            Connection = new TConnectionWithIntervals(Sprinkling, false, true);
            rhoStats->push(Sprinkling->size / Geometry->V);
            double myD1 = (Connection->countD1X() - Connection->countD1Y());
            myStats->push(myD1);
            cout << myStats->getTime();
            delete Connection;
            delete Sprinkling;
        }
        cout << "STATS = " << vecToString(myStats->getStats()) << endl;
        datastream << 2 << "^" << power << "\t";
        datastream << rhoStats->getMean() << "\t";
        datastream << myStats->toString();
        delete Geometry;
    }
}

void testD1_SpatialExponential2_Rectangle() {
    vector<double> parameters(1);
    parameters[0] = 0.5;
    int power = 8;
    //parameters[1] = 0.75 * parameters[0];
    for (int alpha = 1; alpha < 10; alpha++) {
        ofstream datastream("D1-spatialexp-rectangle-dims1x1-alpha-var.txt",
                ios::app);
        int N = pow(2, power);
        double rho = 1.0 / (double) N;
        cout << "N = " << N << endl;
        double cardN1X, cardN1Y, cardD1;
        TGeometry* Geometry;
        TSprinkling* Sprinkling;
        TConnectionWithIntervals* Connection;
        Statistics* myStats = new Statistics("D1");

        Geometry = new TGeometry_SpatialLinear2_Rectangle(1.0, 1.0,
                1.0 / (double) alpha, 1, parameters);

        for (int RUN = 0; RUN < 100; RUN++) {
            Sprinkling = new TSprinkling(rho, Geometry);
            Connection = new TConnectionWithIntervals(Sprinkling, false,
                    true);

            cardN1X = Connection->countD1X();
            cardN1Y = Connection->countD1Y();
            cardD1 = cardN1X - cardN1Y;
            myStats->push(cardD1);
            cout << "N1X, N1Y, D1 = " << cardN1X << ", " << cardN1Y
                    << ", "
                    << cardD1
                    << endl;
            cout << myStats->getTime();
            delete Connection;
            delete Sprinkling;
        }
        delete Geometry;
        cout << "STATS = " << vecToString(myStats->getStats()) << endl;
        datastream << 1 / (double) alpha << "\t";
        datastream << parameters[0] << "\t";
        datastream << 2 << "^" << power << "\t";
        datastream << rho << "\t";
        datastream << myStats->toString();
    }
}

void testD1_SpatialLinear2_Rectangle() {
    vector<double> parameters(1);
    double deltat, deltax;
    double alpha;
    double rho;
    double x0;

    for (int n = 1; n < 2; n++) {
        for (int t = 1; t < 2; t++) {
            for (int a = 1; a < 2; a++) {
                for (int x = 1; x < 4; x++) {
                    for (int r = 0; r < 1; r++) {

                        alpha = 0.9;
                        deltat = 100.0 * alpha;
                        deltax = 10.0 * alpha;
                        x0 = deltax * (double) x / 4.0;
                        parameters[0] = x0;
                        rho = 1.0;

                        ofstream datastream(
                                "D1-spatiallinear-rectangle-100x10-rho1-alpha09-x0var.txt",
                                ios::app);
                        cout << "rho = " << rho << endl;
                        double cardN1X, cardN1Y, cardD1;
                        TGeometry* Geometry;
                        TSprinkling* Sprinkling;
                        TConnectionWithIntervals* Connection;
                        Statistics* myStats = new Statistics("D1", true);

                        Geometry = new TGeometry_SpatialLinear2_Rectangle(
                                (double) deltat, deltax, alpha, 1,
                                parameters);

                        for (int RUN = 0; RUN < 100; RUN++) {
                            Sprinkling = new TSprinkling(Geometry);
                            Sprinkling->writeToFile("mysprinkling.txt");
                            Connection = new TConnectionWithIntervals(
                                    Sprinkling, false, true);

                            cardN1X = Connection->countD1X();
                            cardN1Y = Connection->countD1Y();
                            cardD1 = cardN1X - cardN1Y;
                            myStats->push(cardD1);
                            cout << "N            = " << Sprinkling->size
                                    << endl;
                            cout << "N1X, N1Y, D1 = "
                                    << cardN1X << ", "
                                    << cardN1Y << ", "
                                    << cardD1 << endl;
                            cout << myStats->getTime();
                            delete Connection;
                            delete Sprinkling;
                        }
                        double vol = Geometry->V;
                        delete Geometry;
                        cout << "STATS = "
                                << vecToString(myStats->getStats()) << endl;
                        datastream << rho << "\t";
                        datastream << round(rho * vol) << "\t";
                        datastream << x0 << "\t";
                        datastream << myStats->toString();
                        myStats->dataToFile("mydata.txt");
                    }
                }
            }
        }
    }
}

void testSubregion() {
    vector<double> parameters(1);
    double deltat, deltax;
    double alpha;
    double rho;
    double x0;

    for (int n = 1; n < 2; n++) {
        for (int t = 1; t < 2; t++) {
            for (int a = 1; a < 2; a++) {
                for (int x = 1; x < 2; x++) {
                    for (int r = 0; r < 1; r++) {

                        alpha = 0.9;
                        deltat = 100.0 * alpha;
                        deltax = 10.0 * alpha;
                        x0 = deltax * (double) x / 4.0;
                        parameters[0] = x0;
                        rho = 1.0;

                        ofstream datastream("testSubregion.txt", ios::app);
                        cout << "rho = " << rho << endl;
                        double cardN1X, cardN1Y, cardD1;
                        TGeometry* Geometry;
                        TSprinkling* Sprinkling;
                        TConnectionWithIntervals* Connection;
                        Statistics* myStats = new Statistics("D1", true);

                        Geometry = new TGeometry_SpatialLinear2_Rectangle(
                                (double) deltat, deltax, alpha, 1,
                                parameters);

                        for (int RUN = 0; RUN < 2; RUN++) {
                            Sprinkling = new TSprinkling(Geometry);
                            Sprinkling->writeToFile(
                                    "testSubregion_sprinkling.txt");
                            Connection = new TConnectionWithIntervals(
                                    Sprinkling, false, true);

                            cardN1X = Connection->countD1X();
                            cardN1Y = Connection->countD1Y();
                            cardD1 = cardN1X - cardN1Y;
                            myStats->push(cardD1);
                            cout << "N            = " << Sprinkling->size
                                    << endl;
                            cout << "N1X, N1Y, D1 = "
                                    << cardN1X << ", "
                                    << cardN1Y << ", "
                                    << cardD1 << endl;
                            cout << myStats->getTime();
                            delete Connection;
                            delete Sprinkling;
                        }
                        double vol = Geometry->V;
                        delete Geometry;
                        cout << "STATS = "
                                << vecToString(myStats->getStats()) << endl;
                        datastream << rho << "\t";
                        datastream << round(rho * vol) << "\t";
                        datastream << x0 << "\t";
                        datastream << myStats->toString();
                        myStats->dataToFile("testSubregion_data.txt");
                    }
                }
            }
        }
    }
}

void testD1_Minkowski2_Rectangle_RindlerPartition() {
    vector<double> parameters(2);
    for (int k = 5; k < 6; ++k) {
        parameters[0] = (double) k / 10.0;
        parameters[1] = 0.75 * parameters[0];
        for (int RUNS = 20; RUNS < 30; RUNS++) {
            for (int power = 9; power < 10; power++) {
                ofstream datastream(
                        "D1-rectangle-rindler-dims1x1-power9-k05.txt",
                        ios::app);
                int N = pow(2, power);
                cout << "N = " << N << endl;
                double rho, cardN1X, cardN1Y, cardD1;
                TGeometry* Geometry;
                TSprinkling* Sprinkling;
                TConnectionWithIntervals* Connection;
                Statistics* myStats = new Statistics("D1");

                Geometry = new TGeometry_Minkowski2_Rectangle(1.0, 2.0, 2,
                        parameters);

                for (int RUN = 0; RUN < RUNS * 100; RUN++) {
                    Sprinkling = new TSprinkling(N, false, Geometry);
                    Connection = new TConnectionWithIntervals(Sprinkling,
                            false,
                            true);

                    rho = Sprinkling->size / Geometry->V;
                    cardN1X = Connection->countD1X();
                    cardN1Y = Connection->countD1Y();
                    cardD1 = cardN1X - cardN1Y;
                    myStats->push(cardD1);
//                    cout << "N1X, N1Y, D1 = " << cardN1X << ", " << cardN1Y
//                            << ", "
//                            << cardD1
//                            << endl;
//                    cout << myStats->getTime();
                    delete Connection;
                    delete Sprinkling;
                }
                delete Geometry;
                cout << "STATS = " << vecToString(myStats->getStats())
                        << endl;
                datastream << parameters[0] << "\t";
                datastream << 2 << "^" << power << "\t";
                datastream << rho << "\t";
                datastream << myStats->toString();
            }
        }
    }
}

void testN1_Minkowski2_Rectangle() {
    int N = pow(2, 11);
    vector<double> subradius(1);
    subradius[0];

    TGeometry* Geometry;
    TSprinkling* Sprinkling;
    TConnectionWithIntervals* Connection;
    Statistics* myStats = new Statistics("D1");

    Geometry = new TGeometry_Minkowski2_Rectangle(1, subradius);

    ofstream file("time-comparison.txt");
    file << "d1x\tn1x\n";

    for (int RUN = 0; RUN < 10; RUN++) {
        Sprinkling = new TSprinkling(N, false, Geometry);
        Connection = new TConnectionWithIntervals(Sprinkling, false, true);

        Sprinkling->writeToFileWithSubregion("mysprinkling.txt");
        cout << "subsize = " << Connection->subsize << endl;
        //vector<vector<int> > n1xSets = Connection->getN1X();
        //vector<vector<int> > n1ySets = Connection->getN1XBAR();
        Timer* timer = new Timer("counting D1X");
        cout << "D1X = " << Connection->countD1X() << endl;
        file << timer->toString() << "\t";
//    cout << "N1XX, N1XM = " << Connection->countIntervalsXXOfSize(2) << ", " << Connection->countIntervalsXMOfSize(2) << endl;
        timer = new Timer("counting intervals");
        cout << "N1XX - N1XM = "
                << Connection->countIntervalsXXOfSize(2)
                        - Connection->countIntervalsXMOfSize(2) << endl;
        file << timer->toString() << "\n";
//    for(int i = 0; i < n1xSets.size(); i++) {
//        cout << vecToString(n1xSets[i]) << endl;
//    }
//    cout << "D1Y = " << Connection->countN1XBAR() << endl;
//    cout << "N1YY, N1YM = " << Connection->countIntervalsYYOfSize(2) << ", " << Connection->countIntervalsYMOfSize(2) << endl;
//    cout << "N1YY - N1YM = " << Connection->countIntervalsYYOfSize(2) - Connection->countIntervalsYMOfSize(2) << endl;
//    for(int i = 0; i < n1ySets.size(); i++) {
//        cout << vecToString(n1ySets[i]) << endl;
//    }
        myStats->getTime();
        delete Connection;
        delete Sprinkling;
    }
    delete Geometry;
    cout << "STATS = " << vecToString(myStats->getStats()) << endl;
}

void testOrderIntervals() {
    for (int p = 2; p < 11; p++) {
        int N = pow(2, p);
        int RUNS = 10;
        double subradius = 0.7;
        TGeometry* Geometry = new TGeometry_Minkowski3_Cylinder();

        //ofstream os;
        ofstream os2;
        //os.open("d1cylinder-data.dat", ios::app);
        os2.open("d1cylinder-stat-r=0.3.dat", ios::app);

        double d1mean = 0.0;
        double d1sdev = 0.0;

        for (int RUN = 0; RUN < RUNS; RUN++) {
            cout << "THIS IS RUN NUMBER " << RUN + 1 << endl;
            TSprinkling* Sprinkling = new TSprinkling(N, false, Geometry);
            Sprinkling->writeToFile("testsprinkling140129rsubEQ0.7.txt");
            TConnectionWithIntervals* Connection =
                    new TConnectionWithIntervals(
                            Sprinkling);

            Timer* timer2 = new Timer("countOrderIntervals");
            vector<vector<int> > ints = Connection->getIntervalsOfSize(2);
            timer2->stop();

            int pokingIn = 0;
            int pokingOut = 0;
            for (int i = 0; i < ints.size(); i++) {
                if (!(Sprinkling->inSubregion(ints[i][0]))
                        && !(Sprinkling->inSubregion(ints[i][2]))
                        && (Sprinkling->inSubregion(ints[i][1]))
                        ) pokingIn++;
            }
            for (int i = 0; i < ints.size(); i++) {
                if ((Sprinkling->inSubregion(ints[i][0]))
                        && (Sprinkling->inSubregion(ints[i][2]))
                        && !(Sprinkling->inSubregion(ints[i][1]))
                        ) pokingOut++;
            }
            d1mean += pokingIn - pokingOut;
            d1sdev += pow(pokingIn - pokingOut, 2.0);
            cout << "D1(X,Y) = " << pokingIn - pokingOut << endl;

            delete Connection;
            delete Sprinkling;
        }
        d1mean = d1mean / RUNS;
        d1sdev = sqrt((d1sdev - d1mean * d1mean * RUNS) / (RUNS - 1));
        os2 << subradius << "\t" << RUNS << "\t" << N << "\t"
                << d1mean << "\t"
                << pow(((double) N / Geometry->V), -1 / 3.0) * d1mean
                << "\t"
                << pow(((double) N / Geometry->V), -1 / 3.0) * d1sdev
                << "\n";
        delete Geometry;
    }
}

void testVecIntToString() {
    vector<int> out;
    out.push_back(1);
    out.push_back(2);
    cout << vecToString(out) << endl;
}

void testStatistics() {
    Statistics* myStats = new Statistics("myStats");
    for (int i = 0; i < 50; i++) {
        myStats->push(0.0);
    }
    for (int i = 0; i < 50; i++) {
        myStats->push(10.0);
    }
    cout << myStats->toString();
}

void testDimaDeSitter(double N, double k, int seed) {
    InitRandom();
    double t0 = 1.0;
    TGeometry* Geometry = new TGeometry_dS4Global(t0);
    TSprinkling* Sprinkling = new TSprinkling((int) N, false, Geometry);
    TConnection* Connection = new TConnection(Sprinkling);
    Sprinkling->writeToFile("thissprinkling.dat");
    Connection->writeCausalMatrix("thiscausal.dat");
}

void testConnectionList() {
    TGeometry* Geometry = new TGeometry_Minkowski2_Rectangle();
    TSprinkling* Sprinkling = new TSprinkling(10, Geometry);
    Timer* listTimer = new Timer("list");
    TCauset* Causet = new TCauset(Sprinkling);
    Causet->causalFromSprinkling(Sprinkling);
    cout << Causet->toString();
    listTimer->stop();
    Timer* matrixTimer = new Timer("matrix");
    TConnection* Connection = new TConnection(Sprinkling, false, false);
    cout << Connection->toString();
    matrixTimer->stop();
}

void testPropagators() {
    TGeometry* Geometry = new TGeometry_Minkowski2_Rectangle();
    TSprinkling* Sprinkling = new TSprinkling(10000, Geometry);
    TConnection* Connection = new TConnection(Sprinkling);
    TPropagators* Propagators = new TPropagators(Sprinkling, Connection);
    Timer* timer1 = new Timer("causal");
    Propagators->getCausal();
    timer1->stop();
    Eigen::MatrixXi c = Propagators->getCausal();
    Timer* timer2 = new Timer("retarded");
    Propagators->getRetarded(0.0);
    timer2->stop();
}

//  Newton-Raphson estimation of t0 for given cardinality N and average degree k;
//  NRITERATIONS = 20 gives good results for a range of k/N ~ E-10 to E-3.
double NRt0(double N, double k) {
    double r = k / N;
    double t = 1.3;
    for (int i = 0; i < 20; i++) {
        t = t + (2.0 * (2 + cos(2 * t))
                * (24 * t * pow(cos(t), 4) - cos(t)
                        * (17 + 7 * cos(2 * t)) * sin(t)
                        - 6 * (2 + cos(2 * t)) * log(cos(t)) * sin(2 * t)
                        - 3.0 * PI * r * pow(sin(2 * t) + tan(t), 2)))
                / (3. * (-59 - 34 * cos(2 * t)
                        - 3 * cos(4 * t) + 96 * t * pow(cos(t), 2) / tan(t)
                        - 24 * (2 + cos(2 * t)) * log(cos(t))));
    }
    return t;
}

void debug() {
    TGeometry* Geometry = new TGeometry_EinsteinStaticUniverseD(2);
    TSprinkling* Sprinkling = new TSprinkling(100, false, Geometry);
    delete Sprinkling;
    delete Geometry;
}

int main() {
    InitRandom();
    //testD1_SpatialLinear2_Rectangle();
    //testExtrinsicCurvatureMink3Cylinder();
    //testSpatialLinearSprinkling();
    //testDimaDeSitter(10000, 1, 1);
    //testSpherical();
    //testSprinkling();
    //testConnection();
    //debug();
    //testSprinkling();
    //testExtrinsicCurvatureMink3Cylinder();
    //testExtrinsicCurvatureSmearedMink3();
    testExtrinsicCurvature2MinkowskiDRectangle();
    //testExtrinsicCurvatureTimeSquaredDRectangle();
    cout << "ALL DONE." << endl;
    //cout << NRt0(1.0, 0.0000000001) << endl;
}
