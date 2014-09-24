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
#include "apfMesh.h"
#include "Teuchos_ParameterList.hpp"

/** \namespace awr
  * \brief All AWR error estimation functions
  */
namespace awr {

apf::Field* enrichSolution(apf::Field* sol, const char* name_e);

apf::Field* solveAdjointProblem(
    apf::Mesh* mesh, const Teuchos::ParameterList& params);

}

#endif
