/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVBEZIERSOLUTIONTRANSFER_H
#define CRVBEZIERSOLUTIONTRANSFER_H

/** \file crvBezierShapes.h
  * \brief main file for bezier shape functions */

namespace crv {
class CrvBezierSolutionTransfer
{
  public:
    virtual ~CrvBezierSolutionTransfer();

    virtual bool hasNodesOn(int dimension) = 0;

    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi,
        Entity* vert);

    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities);

    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities);

    virtual void convertToBezierFields(
        int dimension,
        EntityArray &newEntities);
}

CrvBezierSolutionTransfer* createFieldTransfer(f);


#endif
