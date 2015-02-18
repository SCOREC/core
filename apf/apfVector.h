/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFVECTOR_H
#define APFVECTOR_H

/** \file apfVector.h
    \brief The APF linear algebra vector interface */

#include <cmath>
#include "apfArray.h"
#include <iostream>

namespace apf {

/** \brief The mathematical constant pi
  \details although it doesn't fit perfectly in this header,
  there is no more appropriate place for this in APF. */
extern double const pi;

template <std::size_t N>
class Vector;

/** \brief template-generic vector of N doubles
  \details this class implements a linear algebra
  vector whose size is known at compile time.
  It is intended to be used for small, dense vectors
  whose operations benefit from the optimization offered
  by templates and inline methods because they are
  so small that the overhead of function calls and
  conditionals is high.

  For vectors sized at runtime, see apf::DynamicVector.
  For sparse structures or parallel vectors, look
  outside of APF.

  apf::Vector objects should be used in a functional
  programming style, i.e. treated as immutable values
  for the most part, then we expect the compiler to
  optimize the assignments into faster code. */
template <std::size_t N>
class Vector : public Array<double,N>
{
  public:
    /** \brief mandatory */
    Vector() {}
    /** \brief construct from array */
    Vector(double const* v)
    {
      for (std::size_t i=0; i < N; ++i)
        this->elements[i] = v[i];
    }
    /** \brief add two vectors */
    Vector<N> operator+(Vector<N> const& b) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] + b.elements[i];
      return c;
    }
    /** \brief subtract two vectors */
    Vector<N> operator-(Vector<N> const& b) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] - b.elements[i];
      return c;
    }
    /** \brief multiply a vector times a scalar
     \details currently there is no scalar-times-vector operator,
             so do be sure to put scalar on the right hand side. */
    Vector<N> operator*(double s) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] * s;
      return c;
    }
    /** \brief divide a vector by a scalar
      \details equivalent to scaling by 1/s */
    Vector<N> operator/(double s) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] / s;
      return c;
    }
    /** \brief vector dot product
      \details We chose the default vector-vector
      multiplication operator to be the dot product.
      So far this seems to have been a good choice */
    double operator*(Vector<N> const& b) const
    {
      double r=0;
      for (std::size_t i=0; i < N; ++i)
        r += this->elements[i] * b.elements[i];
      return r;
    }
    /** \brief get the vector magnitude */
    double getLength() const {return sqrt((*this)*(*this));}
    /** \brief divide the vector by its magnitude */
    Vector<N> normalize() const {return (*this) / getLength();}
    /** \brief compare two vectors
      \details note that this is an exact comparison,
      so it is not good for things like checking results
      of math routines, which are almost always inexact. */
    bool operator==(Vector<N> const& b) const
    {
      for (std::size_t i=0; i < N; ++i)
        if (this->elements[i] != b.elements[i])
          return false;
      return true;
    }
};

/** \brief 3D vector cross product */
inline Vector<3> cross(Vector<3> const& a, Vector<3> const& b)
{
  Vector<3> r;
  r[0] = a[1]*b[2] - a[2]*b[1];
  r[1] = a[2]*b[0] - a[0]*b[2];
  r[2] = a[0]*b[1] - a[1]*b[0];
  return r;
}

/** \brief vector projection */
template<std::size_t N>
Vector<N> project(Vector<N> const& a, Vector<N> const& b)
{
  return b*((a*b)/(b*b));
}

/** \brief convenience wrapper over apf::Vector<3>
 \details this class adds some functions that could
 not be filled in by templates, mainly
 component-specific initialization and x/y/z names */
class Vector3 : public Vector<3>
{
  public:
    Vector3() {}
    Vector3(Vector<3> const& other):
      Vector<3>(other)
    {}
    /** \brief construct from 3 values
      \details this is commonly used for hardcoding vectors */
    Vector3(double a, double b, double c)
    {
      this->elements[0] = a;
      this->elements[1] = b;
      this->elements[2] = c;
    }
    /** \brief construct from array
     \todo this could maybe be templated */
    Vector3(double const* abc)
    {
      this->elements[0] = abc[0];
      this->elements[1] = abc[1];
      this->elements[2] = abc[2];
    }
    /** \brief write vector to array
     \todo this could be templated */
    void toArray(double* abc) const
    {
      abc[0] = this->elements[0];
      abc[1] = this->elements[1];
      abc[2] = this->elements[2];
    }
    /** \brief read vector from array
     \todo this could be templated */
    void fromArray(const double* abc)
    {
      this->elements[0] = abc[0];
      this->elements[1] = abc[1];
      this->elements[2] = abc[2];
    }
    /** \brief immutable x component */
    double x() const {return this->elements[0];}
    /** \brief immutable y component */
    double y() const {return this->elements[1];}
    /** \brief immutable z component */
    double z() const {return this->elements[2];}
    /** \brief mutable x component */
    double& x() {return this->elements[0];}
    /** \brief mutable y component */
    double& y() {return this->elements[1];}
    /** \brief mutable z component */
    double& z() {return this->elements[2];}
};

}

/** \brief write apf::Vector3 to a C++ stream */
std::ostream& operator<<(std::ostream& s, apf::Vector3 const& v);

#endif
