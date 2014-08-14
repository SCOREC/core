/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFSPR_H
#define APFSPR_H

#include "apf.h"
#include "apfNew.h"
#include "apfDynamicVector.h"
#include "apfDynamicMatrix.h"

namespace apf {

/** @brief compute the gradient of a vector or scalar
  *        field at integration points
  * @param f (In) scalar or vector nodal field
  * @param name (In) name of integration point field
  * @param order (In) integration order of accuracy
  */
Field* getGradIPField(Field* f, const char* name, int order);

/** \brief Computes an isotropic size field from a gradient field.
  *
  * \details This function carries out the majority of the SPR-based
  * error estimator. The input is a gradient-like linear field, which
  * could be strain, stress, or a vector gradient.
  * The output is the linear scalar size field to be used for mesh
  * adaptation.
  *
  * \param eps The gradient-like input field
  * \param adaptRatio The maximum acceptable relative error
  *                   (|err|/|eps|) which determines how fine
  *                   the size field is.
  */
Field* getSPRSizeField(Field* eps, double adaptRatio);

/** \brief Recovers a linear matrix field from an integration point sampling
  */
Field* recoverField(Field* ipMatrixField);

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
