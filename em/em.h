/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef EM_H
#define EM_H


/** \file em.h
 *  \brief The Elegtromagnetics Equilibrated Residual error estimator inteface
 */
#include "apf.h"

#include <mthQR.h>
#include <mth.h>
#include <mth_def.h>

namespace em {


apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g);

apf::Field* equilibrateResiduals(apf::Field* f);










}











#endif
