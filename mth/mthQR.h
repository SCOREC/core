#ifndef MTH_QR_H
#define MTH_QR_H

#include "mthMatrix.h"

/** \file mthQR.h
  * \brief routines for the QR factorization of matrices */

namespace mth {

/** \brief finds the QR decomposition of A
  * \details only 3x3 static and dynamic types
  *          are explicitly instantiated
  * \param a the MxN input matrix (M >= N)
  * \param q the MxM orthogonal output matrix
  * \param r the MxN upper triangular output matrix
  * \returns the rank of A
  */
template <class T, unsigned M, unsigned N>
unsigned decomposeQR(
    Matrix<T,M,N> const& a,
    Matrix<T,M,M>& q,
    Matrix<T,M,N>& r);

/** \brief solves Rx = b for upper triangular R
  * \details when M > N, the lower rows are ignored
  * \param r the MxN (M >= N) upper triangular input matrix
  * \param b the Mx1 right hand side input vector
  * \param x the Nx1 output solution vector
  */
template <class T, unsigned M, unsigned N>
void backsubUT(
    Matrix<T,M,N> const& a,
    Vector<T,M> const& b,
    Vector<T,N>& x);

/** \brief solves Ax = b given A's QR factorization
  * \details when M > N, the least squares problem is solved.
  *          only the dynamic type is explicitly instantiated.
  * \param q the MxM orthogonal input matrix
  * \param r the MxN (M >= N) upper triangular input matrix
  * \param b the Mx1 right hand side input vector
  * \param x the Nx1 output solution vector
  */
template <class T, unsigned M, unsigned N>
void solveFromQR(Matrix<T,M,M> const& q,
    Matrix<T,M,N> const& r,
    Vector<T,M> const& b, Vector<T,N>& x);

/** \brief solves Ax = b using A's QR factorization
  * \details when M > N, the least squares problem is solved.
  *          only the dynamic type is explicitly instantiated.
  * \param a the MxN (M >= N) input matrix
  * \param b the Mx1 right hand side input vector
  * \param x the Nx1 output solution vector
  */
template <class T, unsigned M, unsigned N>
bool solveQR(Matrix<T,M,N> const& a,
    Vector<T,M> const& b, Vector<T,N>& x);

}

#endif
