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
#include "crvTables.h"
#include <cassert>

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

bool repositionEdge(ma::Mesh* m, ma::Entity* tet,
    ma::Entity* edge)
{
  // lets assume we have an edge we want to fix
  // only support second order for now
  int P = m->getShape()->getOrder();
  if (P != 2) return false;

  ma::Entity* verts[4];
  ma::Entity* edges[6];

  ma::Vector pivotPoint;
  ma::Vector edgeVectors[3];
  m->getDownward(tet,0,verts);
  m->getDownward(tet,1,edges);

  // pick a pivotVert, the vertex with the worse jacobian determinant
  ma::Entity* pivotVert;
  int pivotIndex;
  {
    apf::MeshElement* me = apf::createMeshElement(m,tet);

    ma::Entity* edgeVerts[2];
    m->getDownward(edge,0,edgeVerts);
    apf::Matrix3x3 J;
    pivotIndex = apf::findIn(verts,4,edgeVerts[0]);
    assert(pivotIndex >= 0);

    ma::Vector xi = crv::elem_vert_xi[apf::Mesh::TET][pivotIndex];
    apf::getJacobian(me,xi,J);

    double j = apf::getJacobianDeterminant(J,3);
    pivotVert = edgeVerts[0];

    int index = apf::findIn(verts,4,edgeVerts[1]);
    assert(index >= 0);
    xi = crv::elem_vert_xi[apf::Mesh::TET][index];
    apf::getJacobian(me,xi,J);
    if (apf::getJacobianDeterminant(J,3) < j){
      pivotVert = edgeVerts[1];
      pivotIndex = index;
    }
    apf::destroyMeshElement(me);
  }

  m->getPoint(pivotVert,0,pivotPoint);

  // local, of edges around vert, [0,2]
  int edgeIndex = 0;

  for (int i = 0; i < 3; ++i){
    // theres only one point, so reuse this...
    edgeVectors[i] = ma::getPosition(m,edges[vertEdges[pivotIndex][i]])
                   - pivotPoint;
    if (edges[vertEdges[pivotIndex][i]] == edge)
      edgeIndex = i;
  }
  assert(edgeIndex >= 0);
  ma::Vector normal = apf::cross(edgeVectors[(1+edgeIndex) % 3],
      edgeVectors[(2+edgeIndex) % 3]);
  double length = normal.getLength();
  double validity = edgeVectors[edgeIndex]*normal;

  if(validity > 1e-10)
    return false;
  ma::Vector oldPoint = ma::getPosition(m,edge);
  apf::Adjacent adjacent;
  m->getAdjacent(edge,3,adjacent);

  // places the new point at a 20 degree angle with the plane
  ma::Vector newPoint = edgeVectors[edgeIndex] + pivotPoint
      + normal/length*(-validity/length +
          edgeVectors[edgeIndex].getLength()*sin(apf::pi/9.));

  m->setPoint(edge,0,newPoint);

  for (std::size_t i = 0; i < adjacent.getSize(); ++i){
    if (checkValidity(m,adjacent[i],4) > 0){
      m->setPoint(edge,0,oldPoint);
      return false;
    }
  }

  return true;
}

}
