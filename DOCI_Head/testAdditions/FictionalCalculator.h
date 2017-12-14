//
// Created by Wulfix on 04/12/2017.
//

#ifndef DOCIPROJECT_FICTIONALCALCULATOR_H
#define DOCIPROJECT_FICTIONALCALCULATOR_H


#include <abstract/StaticWrapper.h>

class FictionalCalculator : public StaticWrapper {
public:
    double calculateOverlap(unsigned long site1, unsigned long site2) override;

    double calculateOverlap(unsigned long site1, unsigned long site2, unsigned long site3, unsigned long site4) override;
};


#endif //DOCIPROJECT_FICTIONALCALCULATOR_H
