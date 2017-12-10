//
// Created by Wulfix on 02/12/2017.
//

#ifndef DOCIPROJECT_INTEGRALCALCULATOR_H
#define DOCIPROJECT_INTEGRALCALCULATOR_H


class StaticWrapper {
public:
    virtual double calculateOverlap(int site1, int site2)=0;
    virtual double calculateOverlap(int site1, int site2, int site3, int site4)=0;


};


#endif //DOCIPROJECT_INTEGRALCALCULATOR_H
