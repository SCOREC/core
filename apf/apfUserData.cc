/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfUserData.h"

namespace apf {

template <> class UserDataBase<double>;
template <> class UserDataBase<int>;
template <> class UserDataBase<long>;
template <> class UserDataBase<size_t>;

}
