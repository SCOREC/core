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

template <class T, unsigned M, unsigned N>
Matrix<T,M,N> transpose(Matrix<T,M,N> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  Matrix<T,M,N> r(m,n);
  for (unsigned i=0; i < m; ++i)
  for (unsigned j=0; j < n; ++j)
    r(j,i) = a(i,j);
  return r;
}

}

#endif
