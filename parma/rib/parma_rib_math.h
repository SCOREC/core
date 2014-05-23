#ifndef PARMA_RIB_MATH_H
#define PARMA_RIB_MATH_H

#include <apfMatrix.h>

namespace parma {

void getPrincipalEigenvector(apf::Matrix3x3 const& A, apf::Vector3& v);

}

#endif
