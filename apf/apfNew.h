/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFNEW_H
#define APFNEW_H

#include <cstddef>
#include <algorithm>
#include <new>

namespace apf {

/* this wraps a new [] array to at least prevent memory leaks,
   double free() errors, and the like.
   lacking size data, it has a limited but important niche */
template <class T>
class NewArray
{
  public:
    NewArray():elements(0) {}
    NewArray(std::size_t n):elements(0) {allocate(n);}
    ~NewArray() {deallocate();}
    bool allocated() const {return elements;}
    T& operator[](std::size_t i) {return elements[i];}
    T const& operator[](std::size_t i) const {return elements[i];}
    void deallocate() {delete [] elements; elements=0;}
    void allocate(std::size_t n)
    {
      deallocate();
      elements = new T[n];
    }
    void swap(NewArray<T>& other)
    {
      std::swap(elements,other.elements);
    }
  protected:
    T* elements;
  private:
    NewArray(NewArray<T> const& other);
    NewArray<T>& operator=(const NewArray<T>& other);
};

} //namespace apf

#endif
