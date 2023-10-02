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

/** \file crvSnap.h
  * \brief main file for snapping of curved nodes,
  * based off of maSnap.h and maSnap.cc */

namespace crv {

/** \brief gets parametric location on geometry given t in [0,1] on edge */
void transferParametricOnEdgeSplit(
    apf::Mesh* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p);
/** \brief gets parametric location on geometry given barycentric coordinate
     t on tri
    \details uses two linear splits to find location on a triangle, which
    is more stable and handles degeneracy better than using pure barycentrics */
void transferParametricOnTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p);

/** \brief gets parametric location on geometry given t in [0,1] on edge,
    but splitting in geometric space first and then projecting */
void transferParametricOnGeometricEdgeSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p);
/** \brief gets parametric location on geometry given barycentric coordinate
     t on tri, but splitting in geometric space first and then projecting */
void transferParametricOnGeometricTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p);

}

#endif
