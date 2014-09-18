/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_H
#define AWR_H

/** \file awr.h
 *  \brief The AWR error estimation interface
 */

#include "apf.h"

/** \namespace awr
  * \brief All AWR error estimation functions
  */
namespace awr {

/** \brief p-enrichnment of a finite element solution field
  * \param sol finite element solution field
  * \returns an enriched finite element field
  * \details let p be the order of the original solution field.
  *          the enriched solution field is of order p+1
  */
apf::Field* enrichSolution(apf::Field* sol);

}

#endif
