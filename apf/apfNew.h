/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_NEW_H
#define APF_NEW_H

#include <canNewArray.h>

namespace apf {

template <class T>
class NewArray : public can::NewArray<T> {
  public:
    typedef can::NewArray<T> Base;
    NewArray() {}
    NewArray(std::size_t n):Base(n) {}
};

}

#endif
