/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ARRAY_H
#define APF_ARRAY_H

#include <canArray.h>

namespace apf {

template <class T, std::size_t N>
class Array : public can::Array<T, N> {
};

}

#endif
