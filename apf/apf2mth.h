/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_2_MTH_H
#define APF_2_MTH_H

#include <mthMatrix.h>
#include <apfMatrix.h>

namespace apf {

template <std::size_t M>
mth::Vector<double, M> to_mth(apf::Vector<M> const& a)
{
  mth::Vector<double,M> r;
  for (unsigned i = 0; i < M; ++i)
    r(i) = a[i];
  return r;
}

template <std::size_t M, std::size_t N>
mth::Matrix<double, M, N> to_mth(apf::Matrix<M,N> const& a)
{
  mth::Matrix<double,M,N> r;
  for (unsigned i = 0; i < M; ++i)
  for (unsigned j = 0; j < N; ++j)
    r(i,j) = a[i][j];
  return r;
}

}

#endif
