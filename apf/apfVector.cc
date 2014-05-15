/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfVector.h"

namespace apf {

double const pi = 3.14159265358979323846;

}

std::ostream& operator<<(std::ostream& s, apf::Vector3 const& v)
{
  s << '(' << v[0] << ", " << v[1] << ", " << v[2] << ')';
  return s;
}
