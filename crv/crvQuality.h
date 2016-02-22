/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVQUALITY_H
#define CRVQUALITY_H

#include "crv.h"
#include <apf.h>
#include <cassert>

/** \file crvQuality.h
    \brief main file for quality functions defined outside of crvQuality.cc */

namespace crv {

/** \brief subdivide jacobian det using subdivision matrices
    \details see getBezierJacobianDetSubdivisionCoefficients */
void subdivideBezierEntityJacobianDet(int P, int type,
    apf::NewArray<double>& c, apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes);
/** \brief get matrices used for uniform subdivision, 2^dim matrices,
     unrolled into a double */
void getBezierJacobianDetSubdivisionCoefficients(int P, int type,
    apf::NewArray<double>& c);

/** \brief typedef for table of jacobian det subdivision functions */
typedef void (*SubdivisionFunction)(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes);
/** \brief table of jacobian det subdivision functions */
extern const SubdivisionFunction subdivideBezierJacobianDet[apf::Mesh::TYPES];

/** \brief elevate jacobian det to higher order, used in getQuality */
void elevateBezierJacobianDet(int type, int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes);
}

#endif
