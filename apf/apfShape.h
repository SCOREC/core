/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFSHAPE_H
#define APFSHAPE_H

#include "apf.h"
#include "apfNew.h"
#include <sstream>

namespace apf {

class EntityShape
{
  public:
    virtual ~EntityShape();
    virtual void getValues(Vector3 const& xi, NewArray<double>& values) const = 0;
    virtual void getLocalGradients(Vector3 const& xi, NewArray<Vector3>& grads) const = 0;
    virtual int countNodes() const = 0;
};

class FieldShape
{
  public:
    virtual ~FieldShape();
    virtual EntityShape* getEntityShape(int type) = 0;
    virtual bool hasNodesIn(int dimension) = 0;
    virtual int countNodesOn(int type) = 0;
    virtual int getOrder() = 0;
    virtual void getNodeXi(int type, int node, Vector3& xi);
    virtual const char* getName() const = 0;
};

FieldShape* getLagrange(int order);
FieldShape* getConstant(int dimension);
FieldShape* getIPShape(int dimension, int order);
FieldShape* getVoronoiShape(int dimension, int order);

int countElementNodes(FieldShape* s, int type);
}

#endif
