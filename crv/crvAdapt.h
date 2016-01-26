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

namespace crv {

class Adapt : public ma::Adapt
{
  public:
    Adapt(ma::Input* in);
    ma::Tag* validityTag;
};

/* Configures a shape correction input,
 * Since fixing invalid elements is considered
 * Adaptation
 */
ma::Input* configureShapeCorrection(
    ma::Mesh* m, ma::SizeField* f=0,
    ma::SolutionTransfer* s=0);

void adapt(ma::Input* in);

// uses blending to position interior points
// based on edge locations
void repositionInteriorWithBlended(ma::Mesh* m,
    ma::Entity* e);
// based on XJ Luo's thesis
bool repositionEdge(ma::Mesh* m, ma::Entity* tet,
    ma::Entity* edge);

// Split edges marked with ma::SPLIT
void splitEdges(ma::Adapt* a);
// Validity checking is expensive, so tag them with an
// integer to identify if invalid, and where
long markInvalidEntities(Adapt* a);

/* Use a crv version of these
 * because we don't have bitwise operations
 */
int getFlag(Adapt* a, ma::Entity* e);
void setFlag(Adapt* a, ma::Entity* e, int flag);
void clearFlag(Adapt* a, ma::Entity* e);

/* Use an integer to determine the quality tag
 * 0 -> Not checked
 * 1 -> Okay Quality
 * 2-7 -> Vertices are bad
 * 8-13 -> Edges are bad
 * 14-17 -> Face are bad
 * 18 -> Entity itself is bad, this one is the worst
 */
int getQualityTag(ma::Mesh* m, ma::Entity* e,
    ma::Entity* bdry);

/* Take boundary triangles where two edges on the boundary
 * form an angle of 180 (or greater) at a vertex and
 * split the edge opposite them.
 */
int fixLargeBoundaryAngles(Adapt* a);

/* If an edge is flagged as invalid,
 * try and collapse
 * or swap it away
 */
int fixInvalidEdges(Adapt* a);

ma::ShapeHandler* getShapeHandler(ma::Adapt* a);

}

#endif
