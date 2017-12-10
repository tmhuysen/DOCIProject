//
// Created by Wulfix on 04/12/2017.
//


#include "include/DOCI.h"

DOCI::DOCI(unsigned int sites, unsigned int electrons, StaticWrapper& calculator) {



    if (sites <= electrons){
        std::cerr << "Invalid argument: to many electrons";
    }else{
        this->sites = sites;
        this->electrons = electrons;

    }
    auto nbf_ = boost::math::binomial_coefficient<double>(sites, electrons);
    if (nbf_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the sector is too high to be cast into unsigned long.");
    }
    this->nbf = static_cast<unsigned long>(nbf_);
    this->hamiltonian = Eigen::MatrixXd::Zero(this->nbf, this->nbf);

    ad_mat = AddressingMatrix(sites,electrons);
    groundstates = { State {std::numeric_limits<double>::max(),Eigen::VectorXd()} };
    integralCalculator = &calculator;
    calculateDoci(0,1);

    Eigen::EigenSolver<Eigen::MatrixXd> solver(hamiltonian);
    eigenvalues = solver.eigenvalues().real().cast<double>();
    eigenvectors = solver.eigenvectors().real().cast<double>();
    //We extract only the groundstate
    for (int i = 0; i<eigenvalues.size(); i++) {
        groundStates(State {eigenvalues[i], eigenvectors.col(i), 0});
    }

}




void DOCI::calculateDoci(double start, double end) {
    boost::dynamic_bitset<> basic_bit = ad_mat.generateBinaryVector(start * nbf);

    for (unsigned long i = 0; i < nbf * end; i++) {

        for (unsigned long j = 0; j < sites; j++) {


            if (basic_bit[j]){
                addToHamiltonian(integralCalculator->calculateOverlap(j,j),i,i);
                addToHamiltonian(integralCalculator->calculateOverlap(j,j),i,i);

            }
            for(unsigned long l = j; l < sites; l++){

                if(j!=l){
                    boost::dynamic_bitset<> two_target_dia = basic_bit;
                    if (annihilation(two_target_dia, j) && annihilation(two_target_dia, l)){
                        double overlap1 = integralCalculator->calculateOverlap(j,l,j,l);
                        double overlap2 = integralCalculator->calculateOverlap(l,j,l,j);
                        double overlap3 = integralCalculator->calculateOverlap(j,l,l,j);
                        double overlap4 = integralCalculator->calculateOverlap(l,j,j,l);
                        overlap3 *= -1;
                        overlap4 *= -1;
                        addToHamiltonian(overlap1,i,i);
                        addToHamiltonian(overlap2,i,i);
                        addToHamiltonian(overlap3,i,i);
                        addToHamiltonian(overlap4,i,i);

                        addToHamiltonian(overlap1,i,i);
                        addToHamiltonian(overlap2,i,i);
                        addToHamiltonian(overlap3,i,i);
                        addToHamiltonian(overlap4,i,i);


                    }


                }

                boost::dynamic_bitset<> two_target = basic_bit;
                if (annihilation(two_target, j) && creation(two_target, l)){
                    unsigned long address = ad_mat.fetchAddress(two_target);
                    double overlap = integralCalculator->calculateOverlap(j,l,j,l);
                    addToHamiltonian(overlap,i,address);
                    addToHamiltonian(overlap,i,address);
                }





            }



        }

        basic_bit = next_bitset_permutation(basic_bit);

    }
    symmatu(hamiltonian);




}
void DOCI::addToHamiltonian(double value, unsigned long index1, unsigned long index2) {
    hamiltonian(index1, index2) += value;

}










void DOCI::groundStates(State state) {
    if(areSame(state,this->groundstates.at(0))){
        groundstates.push_back(state);
    }
    else{
        if(compareState(state,groundstates.at(0))){
            groundstates = std::vector<State> {state};
        }
    }

}

void DOCI::print() {
    std::cout<<std::endl<<hamiltonian;

}

Eigen::MatrixXd DOCI::getHam() {
    return hamiltonian;
}

const std::vector<State> &DOCI::getGroundstates() const {
    return groundstates;
};





bool compareState(const State &o1, const State &o2) {
    return o1.eigenValue < o2.eigenValue;
}


bool areSame(const State &o1, const State &o2) {
    double precision = 10000000;

    double ELIPSON = (o1.eigenValue > o2.eigenValue) ?  o2.eigenValue/precision : o1.eigenValue/precision   ;

    return fabs(o1.eigenValue - o2.eigenValue) < fabs(ELIPSON);
}
