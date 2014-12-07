/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFDYNAMICVECTOR_H
#define APFDYNAMICVECTOR_H

/** \file apfDynamicVector.h
  \brief Small runtime-sized vectors */

#include "apfDynamicArray.h"
#include "apfVector.h"
#include <math.h>
#include <iostream>

namespace apf {

class DynamicVector;

/** \brief A runtime-sized linear algebra vector of doubles
  \details This class is the runtime-sized equivalent
  of apf::Vector. It keeps all values in a single dynamically
  allocated array. Due to the use of dynamic allocation,
  users should avoid copying this class as much as possible.
  To help with this, we provide things like operator+= instead
  of operator+ to discourage users from creating temporary copies.
  The code for these methods is still inlined in an effort to
  keep your linear algebra running as fast as possible. */

class DynamicVector : public DynamicArray<double>
{
  public:
    /** \brief default constructor - no allocation */
    DynamicVector() {}
    /** \brief construct with (n) allocated elements */
    DynamicVector(std::size_t n):DynamicArray<double>(n) {}
    /** \brief immutable index operator
     \details note that we do inherit the square-bracket
     index operator from DynamicArray, but we use parentheses
     here to be consistent with apf::DynamicMatrix */
    double operator()(std::size_t i) const
    {
      return (*this)[i];
    }
    /** \brief mutable index operator */
    double& operator()(std::size_t i)
    {
      return (*this)[i];
    } 
    /** \brief Add a vector to this vector */
    DynamicVector& operator+=(DynamicVector const& b)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] += b[i];
      return *this;
    }
    /** \brief Subtract a vector from this vector */
    DynamicVector& operator-=(DynamicVector const& b)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] -= b[i];
      return *this;
    }
    /** \brief Multiply this vector by a scalar */
    DynamicVector& operator*=(double s)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] *= s;
      return *this;
    }
    /** \brief Divide this vector by a scalar */
    DynamicVector& operator/=(double s)
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] /= s;
      return *this;
    }
    /** \brief Get the vector dot product */
    double operator*(DynamicVector const& b) const
    {
      double r=0;
      for (std::size_t i=0; i < size; ++i)
        r += (*this)[i] * b[i];
      return r;
    }
    /** \brief Get the vector magnitude */
    double getLength() {return sqrt((*this)*(*this));}
    /** \brief Initialize all elements to zero */
    void zero()
    {
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] = 0;
    }
};

/** \brief convert an apf::Matrix into an apf::DynamicMatrix */
template <std::size_t N>
inline DynamicVector fromVector(Vector<N> other)
{
  DynamicVector result(N);
  for(std::size_t ii = 0; ii < N; ii++)
    result[ii] = other[ii];
  return result;
}

} //namespace apf

std::ostream& operator<<(std::ostream& s, apf::DynamicVector const& x);

#endif
