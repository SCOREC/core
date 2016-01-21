/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "crv.h"
#include "crvAdapt.h"
#include "crvTables.h"
#include <maCoarsen.h>
#include <maEdgeSwap.h>
#include <maOperator.h>
#include <cassert>

/* This is similar to maShape.cc, conceptually, but different enough
 * that some duplicate code makes sense */

namespace crv {

static int vertEdges[4][3] = {{0,2,3},{0,1,4},{1,2,5},{3,4,5}};

static int markEdges(ma::Mesh* m, ma::Entity* e, int tag,
    ma::Entity* edges[3])
{
  int n = 0;
  int md = m->getDimension();
  int type = m->getType(e);
  int vertIndex = tag-2;
  int edgeIndex = tag-8;
  int faceIndex = tag-14;
  int nverts = apf::Mesh::adjacentCount[type][0];
  int nedges = apf::Mesh::adjacentCount[type][1];
  int nfaces = apf::Mesh::adjacentCount[type][2];

  // if we have a single invalid edge, try and swap it
  if(edgeIndex >= 0 && edgeIndex < nedges){
    ma::Downward ed;
    m->getDownward(e,1,ed);
    edges[0] = ed[edgeIndex];
    n = 1;
  } else if (vertIndex >= 0 && vertIndex < nverts){
    //if we have an invalid vertex, try and swap all its edges
    ma::Downward ed;
    m->getDownward(e,1,ed);
    n = md;
    if(md == 2){
      edges[0] = ed[vertIndex];
      edges[1] = ed[(vertIndex+2) % 3];
    } else {
      edges[0] = ed[vertEdges[vertIndex][0]];
      edges[1] = ed[vertEdges[vertIndex][1]];
      edges[2] = ed[vertEdges[vertIndex][2]];
    }
  } else if (faceIndex >= 0 && faceIndex < nfaces){
    ma::Downward ed, faces;
    m->getDownward(e,2,faces);
    m->getDownward(faces[faceIndex],1,ed);
    n = 3;
    edges[0] = ed[0];
    edges[1] = ed[1];
    edges[2] = ed[2];
  }
  return n;
}

class EdgeSwapper : public ma::Operator
{
public:
  EdgeSwapper(Adapt* a)
{
    adapter = a;
    mesh = a->mesh;
    edges[0] = edges[1] = edges[2] = 0;
    simplex = 0;
    edgeSwap = ma::makeEdgeSwap(a);
    md = mesh->getDimension();
    ne = ns = 0;
}
  virtual ~EdgeSwapper()
  {
    delete edgeSwap;
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e)
  {
    int tag = crv::getFlag(adapter,e);

    int ne = markEdges(mesh,e,tag,edges);
    simplex = e;
    return (ne > 0);
  }
  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(edges,1);
  }
  virtual void apply()
  {
    for (int i = 0; i < ne; ++i){
      if (edgeSwap->run(edges[i])){
        ns++;
        crv::clearFlag(adapter,simplex);
      }
    }
  }
private:
  Adapt* adapter;
  ma::Mesh* mesh;
  ma::Entity* simplex;
  ma::Entity* edges[3];
  ma::EdgeSwap* edgeSwap;
  int md;
  int ne;
public:
  int ns;
};

static bool isCornerTriAngleLarge(crv::Adapt *a,
    ma::Entity* tri, int index)
{
  ma::Mesh* m = a->mesh;
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

  // this statement is not exactly a fair comparison, but at least gives
  // some level of control over what is considered an invalid angle
  // an angle is "too large" if the dot product between the corner triangle
  // and the triangle formed by its vertices is negative
  if (cornerNormal*normal < a->input->validQuality)
    return true;

  // one final check to handle the odd case where one of the two tets shared
  // by the angle is invalid at the vertex between the two edges,
  // but all of the faces do not have large angles
  // check both its tets are okay
  if (m->getDimension() == 3 && !isBoundaryEntity(m,tri)){
    ma::Entity* triVerts[3];
    m->getDownward(tri,0,triVerts);
    apf::Up up;
    m->getUp(tri,up);

    apf::Matrix3x3 J;
    for (int i = 0; i < up.n; ++i){
      apf::MeshElement* me = apf::createMeshElement(m,up.e[i]);
      ma::Entity* verts[4];
      m->getDownward(up.e[i],0,verts);

      int tetIndex = apf::findIn(verts,4,triVerts[index]);
      ma::Vector xi = crv::elem_vert_xi[apf::Mesh::TET][tetIndex];
      apf::getJacobian(me,xi,J);
      apf::destroyMeshElement(me);

      if(apf::getJacobianDeterminant(J,3) < a->input->validQuality){
        return true;
      }
    }
  }
  return false;
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
      if(isCornerTriAngleLarge(a,e,(i+1) % 3)){

        ma::Entity* edge = edges[(i+2) % 3];
        if(!ma::getFlag(a,edge,ma::SPLIT)){
          ma::Entity* verts[3];
          m->getDownward(e,0,verts);
          // mark the vertex so linear spaced points are used
          // this is a trick because refine was not made for this
          ma::setFlag(a,verts[(i+1) % 3],ma::BAD_QUALITY);
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

// for each edge, this has the left/right face
static int edgeFaces[6][2] = {{1,0},{2,0},{3,0},{3,1},{1,2},{2,3}};

static ma::Entity* isLargeAngleTet(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;

  ma::Entity* faces[4];
  m->getDownward(e,2,faces);

  int index = -1;

  // find edge that matters
  for (int i = 0; i < 6; ++i){
    if(isBoundaryEntity(m,faces[edgeFaces[i][0]]) &&
        isBoundaryEntity(m,faces[edgeFaces[i][1]])){
      index = i;
      break;
    }
  }
  if(index < 0) return 0;

  apf::FieldShape* fs = m->getShape();
  ma::Entity* edges[6];
  m->getDownward(e,1,edges);

  ma::Entity* edge = 0;

  // lets do a sampling approach. At each point on the edge
  int bt = apf::Mesh::EDGE;
  ma::Entity* leftFace = faces[edgeFaces[index][0]];
  ma::Entity* rightFace = faces[edgeFaces[index][1]];
  apf::MeshElement* leftMe = apf::createMeshElement(m,leftFace);
  apf::MeshElement* rightMe = apf::createMeshElement(m,rightFace);
  apf::Vector3 xi, leftXi,rightXi;

  for (int i = 0; i < fs->countNodesOn(bt); ++i){
    fs->getNodeXi(bt,i,xi);
    leftXi = apf::boundaryToElementXi(m,edges[index],leftFace,xi);
    rightXi = apf::boundaryToElementXi(m,edges[index],rightFace,xi);
    apf::Matrix3x3 leftJ,rightJ;

    apf::getJacobian(leftMe,leftXi,leftJ);
    apf::getJacobian(rightMe,rightXi,rightJ);
    apf::Vector3 leftN = (apf::cross(leftJ[0],leftJ[1])).normalize();
    apf::Vector3 rightN = (apf::cross(rightJ[0],rightJ[1])).normalize();
    // jacobian has two rows, each with a vector
    // these two vectors form the plane
    // compare their directions to see if they are close to coplanar)
    // 10 degree tolerance
    if(std::fabs(leftN*rightN) > 0.9){
      edge = edges[oppEdges[index]];
      break;
    }
  }
  apf::destroyMeshElement(leftMe);
  apf::destroyMeshElement(rightMe);

  return edge;
}

// BAD_QUALITY flag is used on edges to identify them
// as splits for quality, rather than for size refinement
// These two functions handle two seperate situations

// First, triangles are looked at to see if they have an angle > 180 degrees
// There is also a check for the weird situation described above
// Second, tets are looked at to see if they have two faces on the boundary
// where the jacobian determinant is negative along the shared edge,
// indicative of a large angle (curving around a cylinder or sphere).

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
      if (edge)
      {
        assert(m->getType(edge) == 1);
        ma::setFlag(a,edge,ma::SPLIT);
        ma::setFlag(a,edge,ma::BAD_QUALITY);
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
        ma::setFlag(a,edge,ma::BAD_QUALITY);
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
  ma::Entity* edges[3];
  int dimension = m->getDimension();
  do {
    ma::Iterator* it = m->begin(dimension);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      int tag = crv::getFlag(a,e);
      int n = markEdges(m,e,tag,edges);
      for (int i = 0; i < n; ++i){
        ma::Entity* edge = edges[i];
        assert(edge);
        if (edge && !ma::getFlag(a,edge,flag))
        {
          ma::setFlag(a,edge,flag);
          if (a->mesh->isOwned(edge))
            ++count;
        }
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

int fixLargeBoundaryAngles(Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesOppLargeAnglesTet(a);
  count += markEdgesOppLargeAnglesTri(a);

  if ( ! count){
    return 0;
  }
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %d boundary edges with "
      "large angles in %f seconds",count,t1-t0);
  return 0;
}

static void collapseInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesToFix(a,ma::COLLAPSE);
  if ( ! count)
    return;
  ma::Mesh* m = a->mesh;
  int maxDimension = m->getDimension();
  assert(checkFlagConsistency(a,1,ma::COLLAPSE));
  int successCount = 0;
  for (int modelDimension=1; modelDimension <= maxDimension; ++modelDimension)
  {
    checkAllEdgeCollapses(a,modelDimension);
    findIndependentSet(a);
    successCount += ma::collapseAllEdges(a, modelDimension);
  }
  successCount = PCU_Add_Long(successCount);
  double t1 = PCU_Time();
  ma::print("Collapsed %d of %d bad edges "
      "in %f seconds",successCount, count, t1-t0);
}

static void swapInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesToFix(a,ma::BAD_QUALITY);
  if ( ! count)
    return;
  EdgeSwapper es(a);
  ma::applyOperator(a,&es);

  double t1 = PCU_Time();
  ma::print("Swapped %d of %d bad edges "
      "in %f seconds",es.ns, count, t1-t0);
}

int fixInvalidEdges(Adapt* a)
{
  int count = markInvalidEntities(a);
  if ( ! count){
    return 0;
  }
  collapseInvalidEdges(a);
  swapInvalidEdges(a);
  return count;
}

}
