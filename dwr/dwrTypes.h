/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_TYPES_H
#define DWR_TYPES_H

#include <apfArray.h>
#include <Sacado.hpp>

namespace dwr {

typedef Sacado::Fad::DFad<double> AD;
typedef apf::Array<AD,3> AD_Vector3;
typedef apf::Array<AD_Vector3,3> AD_Matrix3x3;

}

#endif
