/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CAN_NEW_ARRAY_H
#define CAN_NEW_ARRAY_H

/** \file canNewArray.h
    \brief wrapper for operator new/delete [] */

#include <cstddef>
#include "canArray.h"

namespace can {

/** \brief wrapper over operator new/delete []
  \details this wrapper is used to automatically
  call delete [] on an array created with new [],
  preventing memory leaks and forming the base
  class for other containers.
  Ideally, all usage of operator new [] should be
  replaced with this. */
template <class T>
class NewArray : public Array<T,0>
{
  public:
    /** \brief default initialize pointer to zero */
    NewArray() {}
    /** \brief construct with (n) elements */
    NewArray(std::size_t n) : Array<T,0>(n) {}
    /** \brief Array destructor frees memory */
    ~NewArray() {}
    /** \brief return true if memory has been allocated */
    bool allocated() const {return this->elems;}
    /** \brief user-callable deallocation helper */
    void deallocate()
    {
      delete [] this->elems;
      this->sz = 0;
      this->elems=0;
    }
    /** \brief user-callable allocation helper
      \details note that no mix of allocate/deallocate
      calls can cause a memory leak */
    void allocate(std::size_t n) {this->resize(n);}
  private:
    NewArray(NewArray<T> const& other);
    NewArray<T>& operator=(const NewArray<T>& other);
};

} //namespace apf

#endif
