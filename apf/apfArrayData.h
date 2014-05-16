/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFARRAYDATA_H
#define APFARRAYDATA_H

#include "apfFieldData.h"

namespace apf {

template <class T>
void freezeFieldData(FieldBase* base);
template <class T>
void unfreezeFieldData(FieldBase* base);
}

#endif

