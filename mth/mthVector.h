/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MTH_VECTOR_H
#define MTH_VECTOR_H

#include "canArray.h"
#include <cmath>
#include <ostream>

/** \file mthVector.h
  * \brief Small compile-time and run-time linear algebra vectors. */

namespace mth {

/** \brief compile-time (static) vector of size N
  * \details This class endows Array<T,N> with the standard
  * mathematical properties of a linear algebra vector. The vector is
  * templated on scalar type so that math can be performed for a variety
  * of (meaningful) scalar types. */
template <class T, unsigned N=0>
class Vector : public can::Array<T,N>
{
  public:
    /** \brief default constructor */
    Vector() {}
    /** \brief construct with n elems
      * \details A dummy constructor Vector(n) is provided so that
      * dynamic and static vectors can be used interchangebly */
    Vector(unsigned n) {}
    /** \brief construct from an array */
    Vector(T const* v)
    {
      for (unsigned i=0; i < N; ++i)
        (*this)[i] = v[i];
    }
    /** \brief mutable index operator
      * \details An index operator (i) is provided so that
      * mth::Vector and mth::Matrix share a common index operator */
    T& operator()(unsigned i) {return (*this[i]);}
    /** \brief immutable index operator */
    T const& operator()(unsigned i) const {return (*this[i]);}
    /** \brief add a vector to this vector */
    Vector<T,N>& operator+=(Vector<T,N> const& b)
    {
      for (unsigned i=0; i < N; ++i)
        (*this)[i] += b[i];
      return (*this);
    }
    /** \brief add two vectors */
    Vector<T,N> operator+(Vector<T,N> const& b) const
    {
      Vector<T,N> r;
      for (unsigned i=0; i < N; ++i)
        r[i] = (*this)[i] + b[i];
      return r;
    }
    /** \brief subtract a vector from this vector */
    Vector<T,N>& operator-=(Vector<T,N> const& b)
    {
      for (unsigned i=0; i < N; ++i)
        (*this)[i] -= b[i];
      return (*this);
    }
    /** \brief subtract two vectors */
    Vector<T,N> operator-(Vector<T,N> const& b) const
    {
      Vector<T,N> r;
      for (unsigned i=0; i < N; ++i)
        r[i] = (*this)[i] - b[i];
      return r;
    }
    /** \brief multiply a vector times a scalar
      * \details currently there is no scalar times vector
      * operator, so do be sure to put scalar on the right
      * hand side */
    Vector<T,N> operator*(T const& s) const
    {
      Vector<T,N> r;
       for (unsigned i=0; i < N; ++i)
         r[i] = (*this)[i] * s;
       return r;
    }
    /** \brief divide a vector by a scalar
      * \details equivalent to scaling by 1/s */
    Vector<T,N> operator/(T const& s) const
    {
      Vector<T,N> r;
       for (unsigned i=0; i < N; ++i)
         r[i] = (*this)[i] / s;
       return r;
    }
    /** \brief vector dot product
      * \details we chose the default vector-vector
      * multiplication operator to be the dot product.
      * so far this seems to have been a good choice */ 
    T operator*(Vector<T,N> const& b) const
    {
      T r = (T)0.0;
      for (unsigned i=0; i < N; ++i)
        r += (*this)[i] * b[i];
      return r;
    }
    /** \brief get the vector magnitude */
    T getLength() const {return sqrt((*this)*(*this));}
    /** \brief divide the vector by its magnitude */
    Vector<T,N> normalize() const {return (*this)/getLength();}
    /** \brief zero the vector */
    void zero()
    {
      for (unsigned i=0; i < N; ++i)
        (*this)[i] = (T)0.0;
    }
};

/** \brief run-time (dynamic) vector
  * \details A runtime-sized equivalent of mth::Vector<T,N>.
  * All values are stored in a single dynamically allocated array.
  * Due to the use of dynamic allocation, users should avoid copying
  * class as much as possible. To help with this, we provide things
  * like operator+= instead of operator+ to discourage users from
  * creating temporary copies. The code for these methods is inlined
  * in an effort to keep linear algebra running quickly */
template <class T>
class Vector<T,0> : public can::Array<T,0>
{
  public:
    /** \brief default constructor - no allocation */
    Vector() {}
    /** \brief construct with n elements */
    Vector(unsigned n) : can::Array<T>(n) {}
    /** \brief mutable index operator
      * \details An index operator (i) is provided so that
      *  mth::Vector and mth::Matrix share a common index operator */
    T& operator()(unsigned i) {return (*this[i]);}
    /** \brief immutable index operator */
    T const& operator()(unsigned i) const {return (*this[i]);}
    /** \brief add a vector to this vector */
    Vector<T,0>& operator+=(Vector<T,0> const& b)
    {
      for (unsigned i=0; i < this->sz; ++i)
        (*this)[i] += b[i];
      return *this;
    }
    /** \brief subtract a vector from this vector */
    Vector<T,0>& operator-=(Vector<T,0> const& b)
    {
      for (unsigned i=0; i < this->sz; ++i)
        (*this)[i] -= b[i];
      return *this;
    }
    /** \brief multiply this vector by a scalar */
    Vector<T,0>& operator*=(T const& s)
    {
      for (unsigned i=0; i < this->sz; ++i)
        (*this)[i] *= s;
      return *this;
    }
    /** \brief divide this vector by a scalar */
    Vector<T,0>& operator/=(T const& s)
    {
      for (unsigned i=0; i < this->sz; ++i)
        (*this)[i] *= s;
      return *this;
    }
    /** \brief get the vector magnitude */
    T getLength() const {return sqrt((*this)*(*this));}
    /** \brief zero the vector */
    void zero()
    {
      for (unsigned i=0; i < this->sz; ++i)
        (*this)[i] = (T)0.0;
    }
};

/** \brief convenience wrapper over apf::Vector<3>
  * \details this class add some functions that could not be
  * filled in by templates, mainly component specific initializaiton
  * and x/y/z names */
template <class T>
class Vector3 : public Vector<T,3>
{
  public:
    /** \brief default constructor */
    Vector3() {}
    /** \brief copy constructor */
    Vector3(Vector<T,3> const& other) : Vector<T,3>(other) {};
    /** \brief construct from 3 values */
    Vector3(T const& a, T const& b, T const& c)
    {
      (*this)[0] = a;
      (*this)[1] = b;
      (*this)[2] = c;
    }
    /** \brief construct from an array */
    Vector3(T const* abc)
    {
      (*this)[0] = abc[0];
      (*this)[1] = abc[1];
      (*this)[2] = abc[2];
    }
    /** \brief write vector to array */
    void toArray(T* abc) const
    {
      abc[0] = (*this)[0];
      abc[1] = (*this)[1];
      abc[2] = (*this)[2];
    }
    /** \brief read vector from array */
    void fromArray(T const* abc)
    {
      (*this)[0] = abc[0];
      (*this)[1] = abc[1];
      (*this)[2] = abc[2];
    }
    /** \brief mutable x component */
    T& x() {return (*this[0]);}
    /** \brief mutable y component */
    T& y() {return (*this[1]);}
    /** \brief mutable z component */
    T& z() {return (*this[2]);}
    /** \brief immutable x component */
    T const& x() const {return (*this)[0];}
    /** \brief immutable y component */
    T const& y() const {return (*this)[1];}
    /** \brief immutable z component */
    T const& z() const {return (*this)[2];}
};

}

template <class T, unsigned N>
std::ostream& operator<<(std::ostream& s, mth::Vector<T,N> const& a)
{
  for (unsigned i=0; i < a.size(); ++i)
    s << a[i] << ' ';
  s << '\n';
  return s;
}

#endif
