/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

/*
 * Here's an experiment. After refinement, re-place interior points based on
 * Boundary entities, rather than the previous structure.
 *
 */
#include "crvAdapt.h"
#include "crvQuality.h"

namespace crv {

void repositionInteriorWithBlended(ma::Mesh* m, ma::Entity* e)
{
  apf::FieldShape * fs = m->getShape();
  int order = fs->getOrder();
  int typeDim = apf::Mesh::typeDimension[m->getType(e)];

  if(!fs->hasNodesIn(typeDim) ||
      getBlendingOrder(apf::Mesh::simplexTypes[typeDim])) return;

  int n = fs->getEntityShape(apf::Mesh::simplexTypes[typeDim])->countNodes();
  int ne = fs->countNodesOn(apf::Mesh::simplexTypes[typeDim]);
  apf::NewArray<double> c;
  getInternalBezierTransformationCoefficients(m,order,1,
      apf::Mesh::simplexTypes[typeDim],c);
  convertInterpolationPoints(m,e,n-ne,ne,c);

}

}
