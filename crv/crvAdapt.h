/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVADAPT_H
#define CRVADAPT_H

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

void snapRefineToBoundary(ma::Adapt* a);
bool repositionInteriorWithBlended(ma::Mesh* m,
    ma::Entity* e);
void repositionInterior(ma::Refine* r);

void splitEdges(ma::Adapt* a);
int markBadQuality(Adapt* a);

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

void fixElementShapes(Adapt* a);

/** \brief experimental function */
ma::ShapeHandler* getShapeHandler(ma::Adapt* a);

void adapt(ma::Input* in);

}

#endif
