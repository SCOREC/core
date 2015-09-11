#ifndef MTH_VECTOR_H
#define MTH_VECTOR_H

#include "mthArray.h"
#include <cmath>

/** \file mthVector.h
  * \brief Small compile-time and run-time linear algebra vectors. */

namespace mth {

/** \brief compile-time (static) vector of size N
  * \details This class endows Array<T,N> with the standard
  * mathematical properties of a linear algebra vector. The vector is
  * templated on scalar type so that math can be performed for a variety
  * of (meaningful) scalar types. */
template <class T, unsigned N=0>
class Vector : public Array<T,N>
{
  public:
    /** \brief default constructor - no allocation */
    Vector() {}
    /** \brief construct from an array */
    Vector(T const* v)
    {
      for (unsigned i=0; i < N; ++i)
        (*this)[i] = v[i];
    }
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
class Vector<T,0> : public Array<T,0>
{
  public:
    /** \brief default constructor - no allocation */
    Vector() {}
    /** \brief construct with n elements */
    Vector(unsigned n) : Array<T>(n) {}
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

}

#endif
