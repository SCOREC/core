/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_DYNAMIC_ARRAY_H
#define APF_DYNAMIC_ARRAY_H

#include <vasDynamicArray.h>

namespace apf {

template <class T>
class DynamicArray : public vas::DynamicArray<T> {
  public:
    typedef vas::DynamicArray<T> Base;
    DynamicArray() {}
    DynamicArray(std::size_t n):Base(n) {}
};

}

#endif
