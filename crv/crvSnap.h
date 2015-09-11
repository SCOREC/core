/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVSNAP_H
#define CRVSNAP_H

#include "apfMesh2.h"
#include "apfShape.h"

/* see maSnap.cc */

namespace crv {

void transferParametricOnEdgeSplit(
    apf::Mesh* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p);
void transferParametricOnTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p);
void transferParametricOnGeometricEdgeSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p);
void transferParametricOnGeometricTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p);

}

#endif
