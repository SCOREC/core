/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVSHAPE_H
#define CRVSHAPE_H

#include "crv.h"
#include "crvAdapt.h"
/** \file crvShape.h
  * \brief main file for shape fixing operations,
  * largely based off of ma functions */

namespace crv {
/** \brief checks if is a boundary entity */
bool isBoundaryEntity(apf::Mesh* m, apf::MeshEntity* e);

/** \brief uses blending to position interior points,
    based on edge locations */
void repositionInteriorWithBlended(ma::Mesh* m,
    ma::Entity* e);

/** \brief Split edges marked with ma::SPLIT and place high order nodes
    using subdivision, see ma::refine */
void splitEdges(ma::Adapt* a);
/** \brief mark invalid entities with validity tag
    \details since validity checking is expensive, do this as little
    as possible and keep information about the check */
int markInvalidEntities(Adapt* a);

/** \brief get validityTag */
int getTag(Adapt* a, ma::Entity* e);
/** \brief set validityTag */
void setTag(Adapt* a, ma::Entity* e, int tag);
/** \brief reset validityTag */
void clearTag(Adapt* a, ma::Entity* e);
/** \brief get validityTag
    \details Use an integer to determine the validity tag
    0 -> Not checked
    1 -> Okay Quality
    2-7 -> Vertices are bad
    8-13 -> Edge is are bad
    14-17 -> Face is are bad
    20 -> Tet itself is bad, this one is the worst

    6*dim + 2 + index */
int getValidityTag(ma::Mesh* m, ma::Entity* e,
    ma::Entity* bdry);

/** \brief Take boundary triangles where two edges on the boundary
    form an angle of 180 (or greater) at a vertex and
    split the edge opposite them. */
int fixLargeBoundaryAngles(Adapt* a);

/** \brief If an edge is flagged as invalid,
    try and collapse or swap it away */
int fixInvalidEdges(Adapt* a);

/** \brief attempts to fix the shape of the
    elements in a same manner as ma::fixElementShape */
void fixCrvElementShapes(Adapt* a);

/** \brief get bezier shape handler */
ma::ShapeHandler* getShapeHandler(ma::Adapt* a);

}

#endif
