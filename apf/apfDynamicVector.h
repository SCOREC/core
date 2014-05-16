/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFDYNAMICVECTOR_H
#define APFDYNAMICVECTOR_H

#include "apfDynamicArray.h"
#include <math.h>

namespace apf {

class DynamicVector;

class DynamicVector : public DynamicArray<double>
{
  public:
    DynamicVector() {}
    DynamicVector(std::size_t n):DynamicArray<double>(n) {}
    double operator()(std::size_t i) const
    {
      return (*this)[i];
    }
    double& operator()(std::size_t i)
    {
      return (*this)[i];
    } 
    DynamicVector& operator+=(DynamicVector const& b)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] += b[i];
      return *this;
    }
    DynamicVector& operator-=(DynamicVector const& b)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] -= b[i];
      return *this;
    }
    DynamicVector& operator*=(double s)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] *= s;
      return *this;
    }
    DynamicVector& operator/=(double s)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] /= s;
      return *this;
    }
    double operator*(DynamicVector const& b) const
    {
      double r=0;
      for (std::size_t i=0; i < size; ++i)
        r += (*this)[i] * b[i];
      return r;
    }
    double getLength() {return sqrt((*this)*(*this));}
    void zero()
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] = 0;
    }
};

} //namespace apf

#endif
