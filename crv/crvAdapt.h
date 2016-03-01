/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVADAPT_H
#define CRVADAPT_H

#include "crv.h"
#include <ma.h>
#include <maAdapt.h>
#include <maInput.h>
#include <maRefine.h>
#include <maTables.h>

/** \file crvAdapt.h
  * \brief main file for curved adaptation, see maAdapt.h */

namespace crv {

/** \brief base crv::Adapt class, looks the same as ma::Adapt,
    but carries tag identifying validity (see crvShape.h) */
class Adapt : public ma::Adapt
{
  public:
    Adapt(ma::Input* in);
    ma::Tag* validityTag;
};

/** \brief change the order of a Bezier Mesh
 * \details going up in order is exact,
 * except for boundary elements, where snapping changes things
 * Going down in order is approximate everywhere
 * */
void changeMeshOrder(apf::Mesh2* m, int newOrder);


}

#endif
