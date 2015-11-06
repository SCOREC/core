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
#include <maShape.h>
#include <cassert>

namespace crv {

static bool repositionInteriorWithBlended(ma::Mesh* m, ma::Entity* e)
{
  apf::FieldShape * fs = m->getShape();

  int type = m->getType(e);
  int P = fs->getOrder();
  int typeDim = apf::Mesh::typeDimension[type];
  int n = fs->getEntityShape(apf::Mesh::simplexTypes[typeDim])->countNodes();
  int ne = fs->countNodesOn(apf::Mesh::simplexTypes[typeDim]);

  assert(typeDim > 1);
  double qualityBefore = getQuality(m,e);

  if(!fs->hasNodesIn(typeDim) ||
      getBlendingOrder(apf::Mesh::simplexTypes[typeDim])) return false;

  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> nodes, newNodes(ne);
  apf::getVectorNodes(elem,nodes);
  apf::destroyElement(elem);

  apf::NewArray<double> c;
  getInternalBezierTransformationCoefficients(m,P,1,
      apf::Mesh::simplexTypes[typeDim],c);

  convertInterpolationPoints(n-ne,ne,nodes,c,newNodes);
  for (int i = 0; i < ne; ++i)
    nodes[n-ne+i] = newNodes[i];

  double qualityAfter = getQuality(type,P,nodes);

  if(qualityAfter < 0.01) return false;
  // should improve things by at least 1% or don't do it.
  if(qualityAfter < qualityBefore*1.01) return false;

  for (int i = 0; i < ne; ++i)
    m->setPoint(e,i,newNodes[i]);
  return true;
}

void repositionInterior(ma::Refine* r)
{
  int successes = 0;
  ma::Mesh* m = r->adapt->mesh;
  for (int d=1; d <= m->getDimension(); ++d){
    for (size_t i=0; i < r->newEntities[d].getSize(); ++i){
      ma::EntityArray& a = r->newEntities[d][i];
      for (size_t i=0; i < a.getSize(); ++i)
        if (apf::getDimension(m,a[i]) == m->getDimension())
          successes += repositionInteriorWithBlended(m,a[i]);
    }
  }
  ma::print("%d successful repositions\n",successes);
}

}
