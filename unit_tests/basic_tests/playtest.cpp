//
// Created by Wulfix on 10/12/2017.
//

#include <testAdditions/FictionalCalculator.h>
#include <include/DOCI.h>
#include "gtest/gtest.h"
#include "include/Extras.h"
#include "include/RHFWrapper.h"
#include <unsupported/Eigen/CXX11/Tensor>

TEST(play_test, test){
    // Specify some data
    const std::string xyzfilename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/DOCIProject/unit_tests/basic_tests/reference/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    double threshold = 1.0e-06;
    std::string basis_name = "STO-3G";

    libwrp::Molecule water (xyzfilename);
    libwrp::Basis basis (water, basis_name);


    hf::rhf::RHF rhf (basis, threshold);
    std::cout<<rhf.C_canonical.innerSize();
    Eigen::MatrixXd check_p = hf::rhf::calculate_P(rhf.C_canonical, 10);
    std::cout<<check_p;
    Eigen::MatrixXd O = Eigen::MatrixXd::Zero (7, 7);
    O.topLeftCorner(5, 5) = 2 * Eigen::MatrixXd::Identity (5, 5);
    std::cout<<std::endl<<O;
    check_p = rhf.C_canonical*O*rhf.C_canonical.adjoint();
    std::cout<<std::endl<<check_p;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl<< " WEW LAD ";
    std::cout<<std::endl<<basis.S;
    std::cout<<std::endl<<basis.V;
    std::cout<<std::endl<<basis.T;
    double sum = 0;
    double el = 0;

    for (int i = 0; i<rhf.C_canonical.innerSize(); i++) {
        for (int j = 0; j<rhf.C_canonical.innerSize(); j++) {
            sum += rhf.C_canonical.col(0)[i]*rhf.C_canonical.col(1)[j]*basis.S(i,j);
            el += rhf.C_canonical.col(0)[i]*rhf.C_canonical.col(1)[j]*basis.T(i,j);


        }

    }
    std::cout<<std::endl<<sum<< "THIS SHOULD BE 0";
    std::cout<<std::endl<<el<< "THIS SHOULD BE NO IDEA";







}