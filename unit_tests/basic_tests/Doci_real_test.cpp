//
// Created by Wulfix on 06/4/2017.
//

#include <testAdditions/FictionalCalculator.h>
#include <include/DOCI.h>
#include "gtest/gtest.h"
#include "include/Extras.h"


TEST(doci_test, test){
    FictionalCalculator fic = FictionalCalculator();
    DOCI testDoci = DOCI(4,2,fic);
    Eigen::MatrixXd testmat(6,6);
    testmat << 4,  2,  2 , 2 , 2 , 0, 2 ,4  ,2 , 2 , 0,  2, 2 , 2 ,4 , 0 , 2,  2, 2 , 2,  0 ,4 , 2,  2, 2 , 0 , 2 , 2 ,4 , 2, 0 , 2 , 2 , 2 , 2 ,4 ;
    testDoci.print();
    Eigen::MatrixXd resultMat = testDoci.getHam();
    ASSERT_TRUE(compareMat(testmat,resultMat));



}