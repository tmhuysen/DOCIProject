//
// Created by Wulfix on 10/12/2017.
//

#ifndef DOCIPROJECT_CONVOLUTEDWRAPPER_H
#define DOCIPROJECT_CONVOLUTEDWRAPPER_H

#include <abstract/StaticWrapper.h>
#include <hf.hpp>
#include <libwrp.hpp>
class ConvolutedWrapper : public StaticWrapper {
private:
    hf::rhf::RHF& rhf_basis;
public:
    explicit ConvolutedWrapper(hf::rhf::RHF& rhf_basis);

    double calculateOverlap(unsigned long site1, unsigned long site2) override;

    double calculateOverlap(unsigned long site1, unsigned long site2, unsigned long site3, unsigned long site4) override;
};



#endif //DOCIPROJECT_RHFWRAPPER_H
