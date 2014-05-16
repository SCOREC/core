/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFDYNAMICARRAY_H
#define APFDYNAMICARRAY_H

#include "apfNew.h"

namespace apf {

/* a dynamically allocated array with size knowledge, allowing
   the full range of copying, resizing, etc.
   This is what std::vector should always have been. */
template <class T>
class DynamicArray
{
  public:
    DynamicArray():size(0) {}
    DynamicArray(std::size_t n):newArray(n),size(n) {}
    DynamicArray(DynamicArray<T> const& other) {copy(other);}
    T& operator[](std::size_t i) {return newArray[i];}
    T const& operator[](std::size_t i) const {return newArray[i];}
    DynamicArray<T>& operator=(DynamicArray<T> const& other)
    {
      if (&other == this) return *this;
      copy(other);
      return *this;
    }
    std::size_t getSize() const {return size;}
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
    void append(T const& v)
    {
      setSize(size+1);
      (*this)[size-1] = v;
    }
    void append(DynamicArray<T> const& other)
    {
      std::size_t oldSize = size;
      setSize(oldSize+other.size);
      for (std::size_t i=0; i < other.size; ++i)
        (*this)[oldSize+i] = other[i];
    }
    typedef T* iterator;
    iterator begin() {return &((*this)[0]);}
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
