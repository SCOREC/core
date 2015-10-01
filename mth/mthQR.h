#ifndef MTH_QR_H
#define MTH_QR_H

#include "mthMatrix.h"

namespace mth {

template <class T, unsigned M, unsigned N>
unsigned decomposeQR(
    Matrix<T,M,N> const& a,
    Matrix<T,M,M>& q,
    Matrix<T,M,N>& r);

}

#endif
