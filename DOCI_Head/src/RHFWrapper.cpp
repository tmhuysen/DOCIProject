//
// Created by Wulfix on 10/12/2017.
//

#include "include/RHFWrapper.h"

RHFWrapper::RHFWrapper(hf::rhf::RHF &rhf_basis):rhf_basis(rhf_basis) {


}

double RHFWrapper::calculateOverlap(int site1, int site2) {
    double kin = 0;
    double nuc = 0;
    for (int i = 0; i<rhf_basis.C_canonical.innerSize(); i++) {
        for (int j = 0; j<rhf_basis.C_canonical.innerSize(); j++) {
            kin += rhf_basis.C_canonical.col(site1)[i]*rhf_basis.C_canonical.col(site2)[j]*rhf_basis.basis.T(i,j);
            nuc += rhf_basis.C_canonical.col(site1)[i]*rhf_basis.C_canonical.col(site2)[j]*rhf_basis.basis.V(i,j);


        }

    }
    return kin+nuc;
}

double RHFWrapper::calculateOverlap(int site1, int site2, int site3, int site4) {
    double rep = 0;
    Eigen::VectorXd so1 = rhf_basis.C_canonical.col(site1);
    Eigen::VectorXd so2 = rhf_basis.C_canonical.col(site2);
    Eigen::VectorXd so3 = rhf_basis.C_canonical.col(site3);
    Eigen::VectorXd so4 = rhf_basis.C_canonical.col(site4);
    for (int i = 0; i<so1.innerSize(); i++) {
        for (int j = 0; j<so1.innerSize(); j++) {
            for (int l = 0; l<so1.innerSize(); l++) {
                for (int k = 0; k<so1.innerSize(); k++) {
                    rep += so1[i]*so2[k]*so3[j]*so4[l]*rhf_basis.basis.tei(i,j,l,k);



                }
            }




        }

    }
    return rep;
}
