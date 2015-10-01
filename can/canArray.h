/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CAN_ARRAY_H
#define CAN_ARRAY_H

/** \file canArray.h
    \brief compile-time size array */

#include <cstddef>

/** \namespace can
  * \brief All CAN symbols are contained in this namespace.
  */
namespace can {

/** \brief compile-time (static) array of size N */
template <class T, unsigned N=0>
class Array
{
  public:
    /** \brief default constructor - necessary */
    Array() {}
    /** \brief copy constructor */
    Array(Array<T,N> const& other) {copy(other);}
    /** \brief elems is destroyed automatically */
    ~Array() {}
    /** \brief assignment operator */
    Array<T,N>& operator=(Array<T,N> const& other)
    {
      copy(other);
      return *this;
    }
    /** \brief mutable index operator */
    T& operator[](unsigned i) {return elems[i];}
    /** \brief immutable index operator */
    T const& operator[](unsigned i) const {return elems[i];}
    /** \brief get the size of this array */
    unsigned size() const {return N;}
    void resize(unsigned n) { (void) n; }
  protected:
    T elems[N];
    void copy(Array<T,N> const& other)
    {
      for (unsigned i=0; i < N; ++i)
        elems[i] = other.elems[i];
    }
};

/** \brief run-time (dynamic) array */ 
template <class T>
class Array<T,0>
{
  public:
    /** \brief default constructor - no allocation */
    Array() : sz(0), elems(0) {}
    /** \brief construct with n elems */
    Array(unsigned n) : sz(0), elems(0) {resize(n);}
    /** \brief copy constructor */
    Array(Array<T,0> const& other) : sz(0), elems(0) {copy(other);}
    /** \brief destructor - need to delete elems */
    ~Array() {delete [] elems;}
    /** \brief assignment operator */
    Array<T,0>& operator=(Array<T,0> const& other)
    {
      copy(other);
      return *this;
    }
    /** \brief mutable index operator */
    T& operator[](unsigned i) {return elems[i];}
    /** \brief immutable index operator */
    T const& operator[](unsigned i) const {return elems[i];}
    /** \brief get the size of this array */
    unsigned size() const {return sz;}
    /** \brief resize the array */
    void resize(unsigned n)
    {
      if (sz == n)
        return;
      delete [] elems;
      sz = n;
      elems = new T[n];
    }
    typedef T* iterator;
    typedef T const* const_iterator;
    T const* begin() const {return elems;}
    T const* end() const {return elems + sz;}
    T* begin() {return elems;}
    T* end() {return elems + sz;}
  protected:
    unsigned sz;
    T* elems;
    void copy(Array<T,0> const& other)
    {
      if (sz != other.sz)
        resize(other.sz);
      for (unsigned i = 0; i < sz; ++i)
        elems[i] = other.elems[i];
    }
};

}

#endif
