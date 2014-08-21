/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFDYNAMICARRAY_H
#define APFDYNAMICARRAY_H

/** \file apfDynamicArray.h
    \brief what most std::vectors should be */

#include "apfNew.h"

namespace apf {

/** \brief a dynamically allocated array with size knowledge.
  \details adding a size variables allowas
   the full range of copying, resizing, etc.
   over the simple apf::NewArray.
   We use composition instead of inheritance
   to prevent exposing things like apf::NewArray::allocate
 */
template <class T>
class DynamicArray
{
  public:
    /** \brief default constructor */
    DynamicArray():size(0) {}
    /** \brief construct with (n) elements */
    DynamicArray(std::size_t n):newArray(n),size(n) {}
    /** \brief copy constructor */
    DynamicArray(DynamicArray<T> const& other) {copy(other);}
    /** \brief mutable index operator */
    T& operator[](std::size_t i) {return newArray[i];}
    /** \brief immutable index operator */
    T const& operator[](std::size_t i) const {return newArray[i];}
    /** \brief assignment operator */
    DynamicArray<T>& operator=(DynamicArray<T> const& other)
    {
      if (&other == this) return *this;
      copy(other);
      return *this;
    }
    /** \brief get size.
        \details there is already a member called "size",
        hence the get/set functions */
    std::size_t getSize() const {return size;}
    /** \brief resize the array */
    void setSize(std::size_t newSize)
    {
      if (size == newSize) return;
      NewArray<T> newArray2(newSize);
      std::size_t commonSize = std::min(size,newSize);
      for (std::size_t i=0; i < commonSize; ++i)
        newArray2[i] = (*this)[i];
      newArray.swap(newArray2);
      size = newSize;
    }
    /** \brief slow element append
        \details this is the one operation where
        std::vector is better. This function is
        here for convenience, but it is O(N) */
    void append(T const& v)
    {
      setSize(size+1);
      (*this)[size-1] = v;
    }
    /** \brief append an array
      \details this is slightly optimized
      for appending the contents of another array. */
    void append(DynamicArray<T> const& other)
    {
      std::size_t oldSize = size;
      setSize(oldSize+other.size);
      for (std::size_t i=0; i < other.size; ++i)
        (*this)[oldSize+i] = other[i];
    }
    /** \brief STL-style iterator type */
    typedef T* iterator;
    /** \brief STL-style begin function */
    iterator begin() {return &((*this)[0]);}
    /** \brief STL-style end function */
    iterator end() {return begin()+size;}
  protected:
    void copy(DynamicArray<T> const& other)
    {
      size = other.size;
      newArray.allocate(size);
      for (std::size_t i=0; i < size; ++i)
        (*this)[i] = other[i];
    }
    NewArray<T> newArray;
    std::size_t size;
};

}//namespace apf

#endif
