/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "crv.h"
#include "crvAdapt.h"
#include <maCoarsen.h>
#include <maEdgeSwap.h>
#include <maOperator.h>
#include <cassert>

/* This is similar to maShape.cc, conceptually, but different enough
 * that some duplicate code makes sense */

namespace crv {

static ma::Entity* isEdgeQualityTag(ma::Mesh* m,
    ma::Entity* e, int tag)
{
  int index = tag-8;
  int ne = apf::Mesh::adjacentCount[m->getType(e)][1];
  if(index >= 0 && index < ne){
    ma::Downward edges;
    m->getDownward(e,1,edges);
    return edges[index];
  }
  return 0;
}

class EdgeSwapper : public ma::Operator
{
  public:
    EdgeSwapper(Adapt* a)
    {
      adapter = a;
      mesh = a->mesh;
      edges[0] = 0;
      edgeSwap = ma::makeEdgeSwap(a);
      md = mesh->getDimension();
      ns = 0;
    }
    virtual ~EdgeSwapper()
    {
      delete edgeSwap;
    }
    virtual int getTargetDimension() {return md;}
    virtual bool shouldApply(ma::Entity* e)
    {
      int tag = crv::getFlag(adapter,e);
      ma::Entity* edge = isEdgeQualityTag(mesh,e,tag);
      if (! edge)
        return false;
      simplex = e;
      edges[0] = edge;
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return o->requestLocality(edges,1);
    }
    virtual void apply()
        {
          if (edgeSwap->run(edges[0])){
            ns++;
            crv::clearFlag(adapter,simplex);
          }
        }
  private:
    Adapt* adapter;
    ma::Mesh* mesh;
    ma::Entity* simplex;
    ma::Entity* edges[1];
    ma::EdgeSwap* edgeSwap;
    int md;
    int ns;
};

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

/* Checks if an angle of a triangle is large (>= 180 degrees)
 * which can be caused by two edges on the boundary curved to it
 *
 */

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
        // mark the vertex so linear spaced points are used
        // this is a trick because refine was not made for this
        ma::setFlag(a,verts[(i+1) % 3],ma::BAD_QUALITY);
        ma::Entity* edge = edges[(i+2) % 3];
        if(!ma::getFlag(a,edge,ma::SPLIT)){
          return edge;
        }
      }
    }
  }
  return 0;
}

/* Checks if an angle of a tet is large (>= 180 degrees)
 * which can be caused by two edges on the boundary curved to it
 *
 * An analytic approach, looking at the control point net of points
 * by comparing surface normals of each pair of control net points
 * adjacent to the triangle is an incredibly complex ordering exercise
 * Rather than attempt to do the ordering, sampling the Jacobian at
 * P+1 points is used. A validity check on this edge could also be used
 */

// we want to split the edge opposite the boundary edge
static int oppEdges[6] = {5,3,4,1,2,0};

static int edgeFaces[6][2] = {{0,1},{0,2},{0,3},{1,3},{1,2},{2,3}};

static ma::Entity* isLargeAngleTet(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;

  ma::Entity* faces[4];
  m->getDownward(e,2,faces);

  int index = -1;

  // find edge that matters
  for (int i = 0; i < 6; ++i)
    if(isBoundaryEntity(m,faces[edgeFaces[i][0]]) &&
        isBoundaryEntity(m,faces[edgeFaces[i][1]])){
      index = i;
      break;
    }

  if(index < 0) return 0;

  apf::MeshElement* me = apf::createMeshElement(m,e);
  apf::FieldShape* fs = m->getShape();
  apf::Matrix3x3 J;
  ma::Entity* edges[6];
  m->getDownward(e,1,edges);
  apf::Vector3 xi, exi;

  ma::Entity* edge = 0;

  int bt = apf::Mesh::EDGE;
  for (int i = 0; i < fs->countNodesOn(bt); ++i){
    fs->getNodeXi(bt,i,xi);
    exi = apf::boundaryToElementXi(m,edges[index],e,xi);
    apf::getJacobian(me,exi,J);
    if(apf::getJacobianDeterminant(J,3) < a->input->validQuality){
      edge = edges[oppEdges[index]];
      ma::Entity* verts[3];
      m->getDownward(edge,0,verts);
      // mark the vertex so linear spaced points are used
      // this is a trick because refine was not made for this
      ma::setFlag(a,verts[0],ma::BAD_QUALITY);
      ma::setFlag(a,verts[1],ma::BAD_QUALITY);
      break;
    }
  }
  apf::destroyMeshElement(me);
  return edge;
}


static int markEdgesOppLargeAnglesTri(Adapt* a)
{
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  do {
    ma::Iterator* it = m->begin(2);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      if(!hasTwoEntitiesOnBoundary(m,e,1)) continue;
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

static int markEdgesOppLargeAnglesTet(Adapt* a)
{
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  do {
    ma::Iterator* it = m->begin(3);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      ma::Entity* edge = isLargeAngleTet(a,e);
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

/* The whole idea is to do the quality check once,
 * and then use the results to mark edges, etc for
 * fixing
 */
static int markEdgesToFix(Adapt* a, int flag)
{

  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  int dimension = m->getDimension();
  do {
    ma::Iterator* it = m->begin(dimension);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      int tag = crv::getFlag(a,e);
      ma::Entity* edge = isEdgeQualityTag(m,e,tag);
      if (edge && !ma::getFlag(a,edge,flag))
      {
        assert(m->getType(edge) == 1);
        ma::setFlag(a,edge,flag);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

void fixLargeBoundaryAngles(Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesOppLargeAnglesTri(a);
  if (a->mesh->getDimension() == 3)
    count += markEdgesOppLargeAnglesTet(a);
  if ( ! count){
    ma::print("no boundary edges with big angles");
    return;
  }
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %d boundary edges with "
      "large angles in %f seconds",count,t1-t0);
}

static void collapseInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  int nedges = markEdgesToFix(a,ma::COLLAPSE);
  if ( ! nedges)
    return;
  ma::Mesh* m = a->mesh;
  int maxDimension = m->getDimension();
  assert(checkFlagConsistency(a,1,ma::COLLAPSE));
  long successCount = 0;
  for (int modelDimension=1; modelDimension <= maxDimension; ++modelDimension)
  {
    checkAllEdgeCollapses(a,modelDimension);
    findIndependentSet(a);
    successCount += ma::collapseAllEdges(a, modelDimension);
  }
  successCount = PCU_Add_Long(successCount);
  double t1 = PCU_Time();
  ma::print("Collapsed %ld bad edges "
      "in %f seconds",successCount, t1-t0);
}

static void swapInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  long count = markBadQuality(a);
  if ( ! count)
    return;
  EdgeSwapper es(a);
  ma::applyOperator(a,&es);

  double t1 = PCU_Time();
  ma::print("Swapped %ld bad edges "
      "in %f seconds",count, t1-t0);
}

void fixInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  long count = markBadQuality(a);
  if ( ! count){
    ma::print("no bad quality elements found");
    return;
  }
  collapseInvalidEdges(a);
  swapInvalidEdges(a);
  double t1 = PCU_Time();
  ma::print("Found %ld bad quality elements"
      " in %f seconds",count, t1-t0);

}

}
