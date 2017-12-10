//
// Created by Wulfix on 06/4/2017.
//

#include <testAdditions/FictionalCalculator.h>
#include <include/DOCI.h>
#include <include/RHFWrapper.h>
#include "gtest/gtest.h"
#include "include/Extras.h"


TEST(doci_test_2, test){

    using namespace std;
    const std::string xyzfilename = "/Users/wulfix/Desktop/Cursussen_Gent/Thesis/DOCIProject/unit_tests/basic_tests/reference/h2o.xyz";  // Specify the relative path to the input .xyz-file (w.r.t. the out-of-source build directory)
    double threshold = 1.0e-06;
    string basis_name = "STO-3G";

    libwrp::Molecule water (xyzfilename);
    libwrp::Basis basis (water, basis_name);


    hf::rhf::RHF rhf (basis, threshold);
    RHFWrapper rhfWrapper = RHFWrapper(rhf);
    DOCI testDoci = DOCI(7,5,rhfWrapper);
    Eigen::MatrixXd resultMat = testDoci.getHam();
    cout<<endl<<endl<<testDoci.getGroundstates().at(0).eigenValue<< " CHECK CHEKC";




}