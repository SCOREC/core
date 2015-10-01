#ifndef MTH_QR_H
#define MTH_QR_H

#include "mthMatrix.h"

namespace mth {

template <class T, unsigned M, unsigned N>
unsigned decomposeQR(
    Matrix<T,M,N> const& a,
    Matrix<T,M,M>& q,
    Matrix<T,M,N>& r);

template <class T, unsigned M, unsigned N>
void backsubUT(
    Matrix<T,M,N> const& a,
    Vector<T,M> const& b,
    Vector<T,N>& x);

}

#endif
