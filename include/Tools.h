#pragma once

#include <sys/_types/_time_t.h>
#include <cmath>
#include "time.h"
#include <new>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <sys/stat.h>

#include "../lib/mersenne/mersenne.cpp"

using namespace std;

const double PI = M_PI;

CRandomMersenne *RanGen;

void InitRandom() {
    RanGen = new CRandomMersenne((int) time(0));
}

void InitRandom(int seed) {
    RanGen = new CRandomMersenne((int) seed);
}

double rnd() {
    return RanGen->Random();
}
// get random gaussian with mean mu and s.d. sigma
double rndGaussian(double mu, double sigma) {
    double phi = 2 * M_PI * rnd();
    double R = sqrt(2 * log(1 / (1 - rnd()))) * sigma;
    return mu + cos(phi) * R;
}
// get random point on the d-sphere of radius r
// note: also works for the d=0-sphere, giving a uniform
// distribution of values equal to +r and -r
vector<double> rndSpherical(int d, double r) {
    vector<double> y(d + 1);
    double norm = 0;
    for (int i = 0; i < d + 1; i++) {
        y[i] = rndGaussian(0.0, 1.0);
        norm += y[i] * y[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < d + 1; i++) {
        y[i] = r * y[i] / norm;
    }
    return y;
}

 // volume of the d-sphere
double volSphere(int dimension, double radius) {
    if(dimension == 0) return 2.00000;
    if(dimension == 1) return 6.28319 * radius;
    if(dimension == 2) return 12.5664 * pow(radius, 2.0);
    if(dimension == 3) return 19.7392 * pow(radius, 3.0);
    if(dimension == 4) return 26.3189 * pow(radius, 4.0);
    else throw std::invalid_argument("I only know volDSphere(d) for d = 0 to d = 4, dawg!" );
 }

//int round(double x) {
//	return (int) (x + 0.5);
//}

double sqr(double x) {
    return x * x;
}

double min(double a, double b) {
    if (a < b)
    return a;
    return b;
}
;

template<typename T> string vecToString(const std::vector<T> &vec) {
    ostringstream stream;
    stream << "(" << vec[0];
    for (int i = 1; i < vec.size(); i++) {
        stream << ", " << vec[i];
    }
    stream << ")";
    return stream.str();
}
;

const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S-", &tstruct);
    return buf;
};

const string currentDate() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d-", &tstruct);
    return buf;
};

// Function: fileExists
/**
    Check if a file exists
@param[in] filename - the name of the file to check

@return    true if the file exists, else false

*/
bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}
