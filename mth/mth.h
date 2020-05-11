/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MTH_H
#define MTH_H

#include "mthTensor.h"

/** \file mth.h
  * \brief templated math functions */

/** \namespace mth
  * \brief All MTH functions are contained in this namespace */
namespace mth {

/** \brief returns vector cross product */
template <class T>
Vector3<T> cross(Vector<T,3> const& a, Vector<T,3> const& b);

/** \brief returns the cross product matrix for the vector */
template <class T>
Matrix3x3<T> cross(Vector<T,3> const& a);

/** \brief returns vector a projected onto vector b */
template <class T, unsigned N>
Vector<T,N> project(Vector<T,N> const& a, Vector<T,N> const& b);

/** \brief vector rejection */
template <class T, unsigned N>
Vector<T,N> reject(Vector<T,N> const& a, Vector<T,N> const& b);

/** \brief transpose of a static matrix */
template <class T, unsigned M, unsigned N>
void transpose(Matrix<T,M,N> const& a, Matrix<T,N,M>& b);

/** \brief determinant of a 2 by 2 matrix */
template <class T>
T determinant(Matrix<T,2,2> const& a);

/** \brief determinant of a 3 by 3 matrix */
template <class T>
T determinant(Matrix<T,3,3> const& a);

/** \brief invert a 2 by 2 matrix */
template <class T>
Matrix<T,2,2> inverse(Matrix<T,2,2> const& a);

/** \brief invert a 3 by 3 matrix */
template <class T>
Matrix<T,3,3> inverse(Matrix<T,3,3> const& a);

/** \brief trace of a tensor */
template <class T>
T trace(Tensor<T> const& a);

/** \brief Frobenius norm of a tensor */
template <class T>
T norm(Tensor<T> const& a);

/** \brief determinant of a tensor */
template <class T>
T determinant(Tensor<T> const& a);

/** \brief transpose of a tensor */
template <class T>
void transpose(Tensor<T> const& a, Tensor<T>& r);

/** \brief inverse of a tensor */
template <class T>
void inverse(Tensor<T> const& a, Tensor<T>& r);

/** \brief identity tensor */
template <class T>
Tensor<T> eye(unsigned d);

}

#endif
