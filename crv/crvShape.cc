/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "crv.h"
#include "crvAdapt.h"
#include <maEdgeSwap.h>
#include <maDoubleSplitCollapse.h>
#include <maShortEdgeRemover.h>
#include <maOperator.h>
#include <cassert>

/* This is similar to maShape.cc, conceptually, but different enough
 * that some duplicate code makes sense */

namespace crv {

static bool isCornerTriAngleLarge(ma::Mesh* m,
    ma::Entity* tri, int index)
{
  apf::Element* elem = apf::createElement(m->getCoordinateField(),tri);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);
  apf::destroyElement(elem);

  ma::Vector normal = ma::getTriNormal(m,tri);

  int P = m->getShape()->getOrder();
  int r = index*(P-1)+3; // index to the right
  int l = ((index+2) % 3)*(P-1)+3+P-2; // index to the left

  ma::Vector cornerNormal = apf::cross((nodes[r]-nodes[index]),
      (nodes[l]-nodes[index]));

  return cornerNormal*normal < 1e-10;

}


static ma::Entity* isLargeAngleTri(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  ma::Entity* edges[3];
  m->getDownward(e,1,edges);
  for (int i = 0; i < 3; ++i)
  {
    ma::Entity* e0 = edges[i];
    ma::Entity* e1 = edges[(i+1) % 3];
    if(isBoundaryEntity(m,e0) && isBoundaryEntity(m,e1))
    {
      if(isCornerTriAngleLarge(m,e,(i+1) % 3)){
        return edges[(i+2) % 3];
      }
    }
  }

  return 0;
}

static ma::Entity* isLargeAngleTet(crv::Adapt* a, ma::Entity* e)

{
  ma::Mesh* m = a->mesh;
  ma::Entity* faces[4];
  m->getDownward(e,2,faces);
  for (int i = 0; i < 4; ++i)
  {
    ma::Entity* edge = isLargeAngleTri(a,faces[i]);
    if(edge) return edge;
  }
  return 0;
};

typedef ma::Entity* (*AngleFunction)(crv::Adapt* a, ma::Entity* e);

const AngleFunction isLargeAngle[4] =
{
    NULL,
    NULL,
    isLargeAngleTri,
    isLargeAngleTet
};

static long markEdgesOppLargeAngles(Adapt* a)
{
  ma::Entity* e;
  long count = 0;
  ma::Mesh* m = a->mesh;
  int dim = m->getDimension();
  ma::Iterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    if(!hasTwoEdgesOnBoundary(m,e)) continue;
    ma::Entity* edge = isLargeAngle[dim](a,e);
    if (edge)
    {
      assert(m->getType(edge) == 1);
      ma::setFlag(a,edge,ma::SPLIT);
      if (a->mesh->isOwned(edge))
        ++count;
    }
  }
  m->end(it);
  return PCU_Add_Long(count);
}


// split large angles at their opposite edge.
static void fixLargeAngles(crv::Adapt* a)
{
  double t0 = PCU_Time();
  long count = crv::markEdgesOppLargeAngles(a);
//  long count = ma::markEdgesToSplit(a);
  if ( ! count) {
    return;
  }
  printf("found %d elements\n",count);
  assert(ma::checkFlagConsistency(a,1,ma::SPLIT));
  ma::Refine* r = a->refine;
  ma::resetCollection(r);
  ma::collectForTransfer(r);
  ma::addAllMarkedEdges(r);
  ma::splitElements(r);
  crv::snapRefineToBoundary(a);
  ma::processNewElements(r);
  ma::destroySplitElements(r);
  crv::repositionInterior(r);
  ma::forgetNewEntities(r);

  double t1 = PCU_Time();
  ma::print("Fixed %li boundary elements in %f seconds",count,t1-t0);
}

void fixElementShapes(crv::Adapt* a)
{
  fixLargeAngles(a);
  return;
}

}
