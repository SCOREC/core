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

namespace apf {

/* \brief Samples the gradient of a scalar or vector field at 
   integration points */
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

/* auxiliary functions for polynomial fitting */
double evalLinearPolynomial(double coefficients[4], Vector3 const& point);
void fitLinearPolynomial(int count,
                         NewArray<Vector3> const& points,
                         NewArray<double> const& values,
                         double coefficients[4]);

Field* recoverGradientByVolume(Field* f);

}

#endif
