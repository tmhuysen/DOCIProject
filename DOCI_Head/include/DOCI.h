//
// Created by Wulfix on 04/12/2017.
//

#ifndef DOCIPROJECT_DOCI_H
#define DOCIPROJECT_DOCI_H


#include <abstract/StaticWrapper.h>
#include <iostream>
#include "include/AddressingMatrix.h"
#include "testAdditions/FictionalCalculator.h"
#include <Eigen/Dense>
#include "include/Extras.h"
struct State {
    double eigenValue;          // The energy of the solution, a.k.a. the eigenvalue
    Eigen::VectorXd eigenVector; // The coefficients of the solution with respect to the given basis, a.k.a. the eigenvector corresponding to the eigenvalue
    int S_z;            // Total projected spin !!!TIMES 2!!!
    //      For total projected spin = 1/2, S_z = 1
    //      For total projected spin = -1, S_z = -2
};

class DOCI {
public:
    DOCI(unsigned int sites, unsigned int electrons, StaticWrapper& calculator);
    Eigen::MatrixXd getHam();
    void print();
private:
    unsigned int sites;
    unsigned int electrons;
    unsigned long nbf;
    Eigen::MatrixXd hamiltonian;
    AddressingMatrix ad_mat;
    StaticWrapper* integralCalculator;
    std::vector<State> groundstates;
public:
    const std::vector<State> &getGroundstates() const;

private:
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    void calculateDoci(double start, double end);

    void addToHamiltonian(double value, unsigned long index1, unsigned long index2);

    void groundStates(State state);
};

bool compareState(const State &o1, const State &o2);
bool areSame(const State &o1, const State &o2);

#endif //DOCIPROJECT_DOCI_H
