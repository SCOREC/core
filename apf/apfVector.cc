/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfVector.h"

namespace apf {

double const pi = 3.14159265358979323846;

}
int factorial(int num)
{
  assert(num < 11);
  static int const table[11] =
    {1,1,2,6,24,120,720,5040,40320,362880,3628800};
  return table[num];
}
int binomial(int n, int i)
{
  static int const table[28] =
  {1,1,1,1,1,1,1,1,2,3,4,5,6,1,3,6,10,15,1,4,10,20,1,5,15,1,6,1};
  return table[i*7 - (i-1)*i/2 + n-i];
}
std::ostream& operator<<(std::ostream& s, apf::Vector3 const& v)
{
  s << '(' << v[0] << ", " << v[1] << ", " << v[2] << ')';
  return s;
}
