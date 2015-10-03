/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_2_MTH_H
#define APF_2_MTH_H

#include <mthVector.h>
#include <apfVector.h>

namespace apf {

template <std::size_t M>
mth::Vector<double, (unsigned)M > to_mth(apf::Vector<M> const& a)
{
  mth::Vector<double,M> r;
  for (unsigned i = 0; i < M; ++i)
    r(i) = a[i];
  return r;
}

}

#endif
