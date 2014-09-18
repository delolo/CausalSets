/*
 * UnitPropagators.h
 *
 *  Created on: 22 Jun 2014
 *      Author: mb909
 */
#pragma once

using namespace std;

#include "UnitSprinkling.h"
#include "UnitConnection.h"
#include "/usr/local/include/eigen3/Eigen/Core"
#include "/usr/local/include/eigen3/Eigen/Dense"

class TPropagators {
private:
    int size;
    double rho;
    Eigen::MatrixXi causal;
    Eigen::MatrixXd retarded;
    Eigen::MatrixXd wightman;
    TConnection* Connection;
    TSprinkling* Sprinkling;

public:
    TPropagators(TSprinkling* Sprinkling, TConnection* Connection) {
        this->Connection = Connection;
        this->size = Connection->size;
        this->Sprinkling = Sprinkling;
        this->rho = this->size / Sprinkling->geometry->V;
        setCausal(Connection);
    }

    Eigen::MatrixXi getCausal() {
        return causal;
    }

    Eigen::MatrixXd getRetarded(double m) {
        if(retarded.size() != 0) return retarded;
        double a = 0.5;
        double b = 0.5 * m * m / rho;
        Eigen::MatrixXi id = Eigen::MatrixXi::Identity(size,size);
        Eigen::MatrixXd A = id.cast<double>() + causal.cast<double>() * a * b;
        Eigen::MatrixXd B = causal.cast<double>() * a;
        this->retarded = A.triangularView<Eigen::Upper>().solve(B);
        return this->retarded;
    }

    Eigen::MatrixXd getWightman() {
        Eigen::MatrixXd delta = retarded - retarded.transpose();
        return wightman;
    }

private:
    void setCausal(TConnection* Connection) {
        causal = Eigen::MatrixXi(size, size);
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < size - 1 - i; j++) {
                causal(i,j) = (int) Connection->causal(i,j);
            }
        }
    }
};
