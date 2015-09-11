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

void subdivideBezierEdgeJacobianDet(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> (&subNodes)[2]);

void subdivideBezierTriangleJacobianDet(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> (&subNodes)[4]);

}

#endif
