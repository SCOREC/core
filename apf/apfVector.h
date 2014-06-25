/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFVECTOR_H
#define APFVECTOR_H

#include <cmath>
#include "apfArray.h"
#include <iostream>

namespace apf {

//not sure where else to put this...
extern double const pi;

template <std::size_t N>
class Vector;

template <std::size_t N>
class Vector : public Array<double,N>
{
  public:
    Vector<N> operator+(Vector<N> const& b) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] + b.elements[i];
      return c;
    }
    Vector<N> operator-(Vector<N> const& b) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] - b.elements[i];
      return c;
    }
    Vector<N> operator*(double s) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] * s;
      return c;
    }
    Vector<N> operator/(double s) const
    {
      Vector<N> c;
      for (std::size_t i=0; i < N; ++i)
        c.elements[i] = this->elements[i] / s;
      return c;
    }
    double operator*(Vector<N> const& b) const
    {
      double r=0;
      for (std::size_t i=0; i < N; ++i)
        r += this->elements[i] * b.elements[i];
      return r;
    }
    double getLength() const {return sqrt((*this)*(*this));}
    Vector<N> normalize() const {return (*this) / getLength();}
    bool operator==(Vector<N> const& b) const
    {
      for (std::size_t i=0; i < N; ++i)
        if (this->elements[i] != b.elements[i])
          return false;
      return true;
    }
};

inline Vector<3> cross(Vector<3> const& a, Vector<3> const& b)
{
  Vector<3> r;
  r[0] = a[1]*b[2] - a[2]*b[1];
  r[1] = a[2]*b[0] - a[0]*b[2];
  r[2] = a[0]*b[1] - a[1]*b[0];
  return r;
}

template<std::size_t N>
Vector<N> project(Vector<N> const& a, Vector<N> const& b)
{
  return b*((a*b)/(b*b));
}

class Vector3 : public Vector<3>
{
  public:
    Vector3() {}
    Vector3(Vector<3> const& other):
      Vector<3>(other)
    {}
    Vector3(double a, double b, double c)
    {
      this->elements[0] = a;
      this->elements[1] = b;
      this->elements[2] = c;
    }
    Vector3(double* abc)
    {
      this->elements[0] = abc[0];
      this->elements[1] = abc[1];
      this->elements[2] = abc[2];
    }
    void toArray(double* abc) const
    {
      abc[0] = this->elements[0];
      abc[1] = this->elements[1];
      abc[2] = this->elements[2];
    }
    void fromArray(const double* abc)
    {
      this->elements[0] = abc[0];
      this->elements[1] = abc[1];
      this->elements[2] = abc[2];
    }
    double x() const {return this->elements[0];}
    double y() const {return this->elements[1];}
    double z() const {return this->elements[2];}
    double& x() {return this->elements[0];}
    double& y() {return this->elements[1];}
    double& z() {return this->elements[2];}
};

}

std::ostream& operator<<(std::ostream& s, apf::Vector3 const& v);

#endif
