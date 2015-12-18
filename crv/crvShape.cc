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
        ma::Entity* verts[3];
        m->getDownward(e,0,verts);
        // mark the vertex
        ma::setFlag(a,verts[(i+1) % 3],ma::SNAP);
        ma::Entity* edge = edges[(i+2) % 3];
        if(!ma::getFlag(a,edge,ma::SPLIT)){
          return edge;
        }
      }
    }
  }

  return 0;
}

static int markEdgesOppLargeAngles(Adapt* a)
{
  ma::Entity* e;
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  do {
    ma::Iterator* it = m->begin(2);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      if(!hasTwoEdgesOnBoundary(m,e)) continue;
      ma::Entity* edge = isLargeAngleTri(a,e);
      if (edge && !ma::getFlag(a,edge,ma::SPLIT))
      {
        assert(m->getType(edge) == 1);
        ma::setFlag(a,edge,ma::SPLIT);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

void fixElementShapes(crv::Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesOppLargeAngles(a);
  if ( ! count)
    return;
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %d boundary edges with "
      "large angles in %f seconds\n",count,t1-t0);
}

}
