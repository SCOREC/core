/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFARRAY_H
#define APFARRAY_H

#include <cstddef>
#include <cassert>

namespace apf {

template <class T, std::size_t N>
class Array
{
  public:
    enum { size = N };
    Array() {}
    Array(Array<T,N> const& other) {copy(other);}
    ~Array() {}
    Array<T,N>& operator=(Array<T,N> const& other)
    {
      copy(other);
      return *this;
    }
    T& operator[](std::size_t i) {return elements[i];}
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
