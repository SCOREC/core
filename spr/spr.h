/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef SPR_H
#define SPR_H

/** \file spr.h */

#include "apf.h"
#include "apfNew.h"
#include "apfDynamicVector.h"
#include "apfDynamicMatrix.h"

namespace spr {

/** @brief compute the gradient of a vector or scalar
  *        field at integration points
  * @param f (In) scalar or vector nodal field
  * @param name (In) name of integration point field
  * @param order (In) integration order of accuracy
  */
apf::Field* getGradIPField(apf::Field* f,
                           const char* name, 
                           int order);

/** @brief recover a nodal field using patch recovery
  * @param f (In) integration point field
  */
apf::Field* recoverField(apf::Field* ip_field);

apf::Field* getSPRSizeField(apf::Field* f, double adapt_ratio);

/** @brief solve linear least squares problem Ax=b
  * @param A (In) mxn matrix
  * @param x (Out) nx1 solution vector
  * @param b (In) mx1 right hand side vector
  */
void solveSVD(apf::DynamicMatrix& A,
              apf::DynamicVector& x,
              apf::DynamicVector& b);

}

#endif
