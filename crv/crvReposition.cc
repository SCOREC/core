/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvBezier.h"
#include "crvQuality.h"
#include "crvShape.h"
#include "crvTables.h"
#include <pcu_util.h>

namespace crv {

void repositionInteriorWithBlended(ma::Mesh* m, ma::Entity* e)
{
  apf::FieldShape * fs = m->getShape();
  int order = fs->getOrder();
  int typeDim = apf::Mesh::typeDimension[m->getType(e)];
  int type = apf::Mesh::simplexTypes[typeDim];

  if(!fs->hasNodesIn(typeDim) || getBlendingOrder(type))
    return;

  int n = fs->getEntityShape(type)->countNodes();
  int ne = fs->countNodesOn(type);
  apf::NewArray<double> c;
  getInternalBezierTransformationCoefficients(m,order,1,type,c);
  convertInterpolationPoints(m,e,n-ne,ne,c);

}

}
