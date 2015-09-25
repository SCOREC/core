/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MTH_H
#define MTH_H

#include "mthMatrix.h"

/** \file mth.h
  * \brief templated math functinos */

/** \namespace mth
  * \brief All MTH functions are contained in this header */
namespace mth {

/** \brief returns vector cross product */
template <class T>
Vector3<T> cross(Vector3<T> const& a, Vector3<T> const& b);

/** \brief returns vector a projected onto vector b */
template <class T, unsigned N>
Vector<T,N> project(Vector<T,N> const& a, Vector<T,N> const& b);

/** \brief vector rejection */
template <class T, unsigned N>
Vector<T,N> reject(Vector<T,N> const& a, Vector<T,N> const& b);

/** \brief return an m by m static or dynamic identity matrix */
template <class T, unsigned M>
Matrix<T,M,M> eye(unsigned m);

/** \brief transpose of a static or dynamic matrix */
template <class T, unsigned M, unsigned N>
Matrix<T,M,N> transpose(Matrix<T,M,N> const& a);

/** \brief trace of a square static or dynamic matrix */
template <class T, unsigned M>
T trace(Matrix<T,M,M> const& a);

/** \brief determinant of a square static or dynamic matrix
  * \details Only 2x2 and 3x3 matrices currently supported. */
template <class T, unsigned M>
T determinant(Matrix<T,M,M> const& a);

/** \brief invert a static or dynamic square matrix
  * \details Only 2x2 currently supported */
template <class T, unsigned M>
Matrix<T,M,M> inverse(Matrix<T,M,M> const& a);

/** \brief get the deviatoric part of a square static or dynamic matrix */
template <class T, unsigned M>
Matrix<T,M,M> deviatoric(Matrix<T,M,M> const& a);

/** \brief get the Frobenius norm of a static or dynamic matrix */
template <class T, unsigned M, unsigned N>
T norm(Matrix<T,N,N> const& a);

}

#include "mth_def.h"

#endif
