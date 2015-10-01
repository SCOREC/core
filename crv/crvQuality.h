/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVQUALITY_H
#define CRVQUALITY_H

#include "crv.h"

namespace crv {

typedef void (*ElevateFunction)(int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes);

typedef void (*SubdivisionFunction)(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes);

extern ElevateFunction elevateBezierJacobianDet[apf::Mesh::TYPES];

extern SubdivisionFunction subdivideBezierJacobianDet[apf::Mesh::TYPES];

}

#endif
