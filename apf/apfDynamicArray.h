/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_DYNAMIC_ARRAY_H
#define APF_DYNAMIC_ARRAY_H

#include <canArray.h>

namespace apf {

template <class T>
class DynamicArray : public can::Array<T, 0> {
  public:
    typedef can::Array<T, 0> Base;
    DynamicArray() {}
    DynamicArray(std::size_t n):Base((unsigned)n) {}
    /** \brief get size.
        \details this is here for backwards compatibility */
    std::size_t getSize() const {return Base::size();}
    /** \brief resize the array
        \details this is here for backwards compatibility
        with apf::DynamicArray, hence the different code
        because it preserves common elements */
    void setSize(unsigned n)
    {
      if (Base::size() == n)
        return;
      T* newElems = new T[n];
      unsigned commonSize = Base::size();
      if (n < Base::size())
        commonSize = n;
      for (unsigned i = 0; i < commonSize; ++i)
        newElems[i] = (*this)[i];
      delete [] Base::elems;
      Base::elems = newElems;
      Base::sz = n;
    }
    /** \brief slow element append
        \details this is the one operation where
        std::vector is better. This function is
        here for convenience, but it is O(N) */
    void append(T const& v)
    {
      setSize(Base::size() + 1);
      (*this)[Base::size() - 1] = v;
    }
    /** \brief append an array
      \details this is slightly optimized
      for appending the contents of another array. */
    void append(DynamicArray<T> const& other)
    {
      std::size_t oldSize = Base::size();
      setSize(oldSize + other.size());
      for (std::size_t i = 0; i < other.size(); ++i)
        (*this)[oldSize + i] = other[i];
    }
};

}

#endif
