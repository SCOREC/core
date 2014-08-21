/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFARRAY_H
#define APFARRAY_H

/** \file apfArray.h
    \brief APF compile-time size array */

#include <cstddef>
#include <cassert>

namespace apf {

/** \brief an array of N items of type T */
template <class T, std::size_t N>
class Array
{
  public:
    /** \brief type::size for convenience */
    enum { size = N };
    /** \brief constructor: no default values */
    Array() {}
    /** \brief copy constructor */
    Array(Array<T,N> const& other) {copy(other);}
    /** \brief element destruction already automatic */
    ~Array() {}
    /** \brief copy array contents.
      \details this is one of the advantages
      of apf::Array over C arrays: they are
      copy-constructible and assignable */
    Array<T,N>& operator=(Array<T,N> const& other)
    {
      copy(other);
      return *this;
    }
    /** \brief mutable index operator */
    T& operator[](std::size_t i) {return elements[i];}
    /** \brief immutable index operator */
    T const& operator[](std::size_t i) const {return elements[i];}
  protected:
    void copy(Array<T,N> const& other)
    {
      for (std::size_t i=0; i < N; ++i)
        elements[i] = other.elements[i];
    }
    T elements[N];
};

}

#endif
