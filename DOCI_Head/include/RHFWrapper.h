//
// Created by Wulfix on 10/12/2017.
//

#ifndef DOCIPROJECT_RHFWRAPPER_H
#define DOCIPROJECT_RHFWRAPPER_H

#include <abstract/StaticWrapper.h>
#include <hf.hpp>
#include <libwrp.hpp>
class RHFWrapper : public StaticWrapper {
private:
    hf::rhf::RHF& rhf_basis;
public:
    explicit RHFWrapper(hf::rhf::RHF& rhf_basis);

    double calculateOverlap(int site1, int site2) override;

    double calculateOverlap(int site1, int site2, int site3, int site4) override;
};



#endif //DOCIPROJECT_RHFWRAPPER_H
