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

void subdivideBezierEntityJacobianDet(int P, int type,
    apf::NewArray<double>& c, apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes);

void getBezierJacobianDetSubdivisionCoefficients(apf::Mesh* m, int P, int type,
    apf::NewArray<double> & c);

typedef void (*ElevateFunction)(int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes);

typedef void (*SubdivisionFunction)(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes);

extern const ElevateFunction elevateBezierJacobianDet[apf::Mesh::TYPES];

extern const SubdivisionFunction subdivideBezierJacobianDet[apf::Mesh::TYPES];

}

#endif
