/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvAdapt.h"
#include "crv.h"
#include <maSnap.h>
#include <PCU.h>
#include <cassert>

namespace crv {

static void snapToBoundary(ma::Adapt* a){

  if ( ! a->input->shouldSnap)
    return;

  ma::Mesh* m = a->mesh;
  ma::Refine* r = a->refine;
  int P = m->getShape()->getOrder();
  apf::FieldShape* fs = m->getShape();
  BezierCurver bc(m,P,0);
  for (int d=1; d <= m->getDimension(); ++d){
    for (size_t i=0; i < r->newEntities[d].getSize(); ++i){
      ma::EntityArray& a = r->newEntities[d][i];
      for (size_t i=0; i < a.getSize(); ++i)

        if(m->getModelType(m->toModel(a[i])) < m->getDimension()){
          snapToInterpolate(m,a[i]);
        }
    }
  }
  for (int d = m->getDimension(); d >=1; --d){
    for (size_t i=0; i < r->newEntities[d].getSize(); ++i){
      ma::EntityArray& a = r->newEntities[d][i];
      for (size_t i=0; i < a.getSize(); ++i)
        if (m->getType(a[i]) != apf::Mesh::VERTEX &&
            m->getModelType(m->toModel(a[i])) < m->getDimension()){
          int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
          int ni = fs->countNodesOn(d);
          apf::NewArray<double> c;
          crv::getBezierTransformationCoefficients(P,d,c);
          bc.convertInterpolationPoints(a[i],n,ni,c);
        }
    }
  }
}

bool refine(ma::Adapt* a)
{
  double t0 = PCU_Time();
  --(a->refinesLeft);
  long count = ma::markEdgesToSplit(a);
  if ( ! count) {
    return false;
  }
  assert(ma::checkFlagConsistency(a,1,ma::SPLIT));
  ma::Refine* r = a->refine;
  ma::resetCollection(r);
  ma::collectForTransfer(r);
  ma::collectForMatching(r);
  ma::addAllMarkedEdges(r);
  ma::splitElements(r);
  crv::snapToBoundary(a);

  ma::processNewElements(r);

  ma::destroySplitElements(r);

  crv::repositionInterior(r);
  ma::forgetNewEntities(r);
  double t1 = PCU_Time();
  ma::print("refined %li edges in %f seconds",count,t1-t0);
  return true;
}

}
