//
// Created by Wulfix on 04/12/2017.
//

#ifndef DOCIPROJECT_FICTIONALCALCULATOR_H
#define DOCIPROJECT_FICTIONALCALCULATOR_H


#include <abstract/StaticWrapper.h>

class FictionalCalculator : public StaticWrapper {
public:
    double calculateOverlap(int site1, int site2) override;

    double calculateOverlap(int site1, int site2, int site3, int site4) override;
};


#endif //DOCIPROJECT_FICTIONALCALCULATOR_H
