/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVSHAPE_H
#define CRVSHAPE_H

#include "crv.h"

namespace crv {
/** \brief checks if is a boundary entity */
bool isBoundaryEntity(apf::Mesh* m, apf::MeshEntity* e);

}

#endif
