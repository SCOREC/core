/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_H
#define DWR_H

/** \page dwr Dual Weighted Residuals
  * This is the main API for SCOREC's DWR library
  *
  * These functions provide capabilities to solve adjoint
  * boundary value problems using a specified functional
  * quantity of interest.
  */

/** \file dwr.h */

#include "dwrElasticityProblem.h"

namespace apf {
class Mesh;
class Field;
}

/** \namespace dwr
  * \brief All DWR symbols */
namespace dwr {

/** \brief create a new elasticity problem
  * \details see dwrElasticityProblem.h */
ElasticityProblem* createElasticityProblem();

/** \brief destroy an elasticity problem
  * \details users should call this once solve has
  * been run for the elasticity problem */
void destroyElasticityProblem(ElasticityProblem* p);

/** \brief estimate error
  * \brief details estimate the error for an
  * elasticity problem
  * \param p elasticity problem object */
void estimateError(ElasticityProblem* p);

}

#endif
