/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFNEW_H
#define APFNEW_H

/** \file apfNew.h
    \brief wrapper for operator new/delete [] */

#include <cstddef>
#include <algorithm>
#include <new>

namespace apf {

/** \brief wrapper over operator new/delete []
  \details this wrapper is used to automatically
  call delete [] on an array created with new [],
  preventing memory leaks and forming the base
  class for other containers.
  Ideally, all usage of operator new [] should be
  replaced with this.
  However, since this class does not store its
  own size, it cannot be copied and its niche
  is limited. see apf::DynamicArray for the
  next step */
template <class T>
class NewArray
{
  public:
    /** \brief default initialize pointer to zero */
    NewArray():elements(0) {}
    /** \brief construct with (n) elements */
    NewArray(std::size_t n):elements(0) {allocate(n);}
    /** \brief destructor automatically frees memory */
    ~NewArray() {deallocate();}
    /** \brief return true if memory has been allocated */
    bool allocated() const {return elements;}
    /** \brief mutable index operator */
    T& operator[](std::size_t i) {return elements[i];}
    /** \brief immutable index operator */
    T const& operator[](std::size_t i) const {return elements[i];}
    /** \brief user-callable deallocation helper
        \details remember that operator delete [] is
                 a no-op on a zero pointer */
    void deallocate() {delete [] elements; elements=0;}
    /** \brief user-callable allocation helper
      \details note that no mix of allocate/deallocate
      calls can cause a memory leak */
    void allocate(std::size_t n)
    {
      deallocate();
      elements = new T[n];
    }
    /** \brief swap pointers with another NewArray */
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
