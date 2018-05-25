#include "apf.h"
#include "apfMesh.h"
#include "apfShape.h"
#include <iostream>
#include <pcu_util.h>

namespace apf {

static Matrix3x3 getJacobianInMetric(const Matrix3x3& J, const Matrix3x3& Q)
{
  Matrix3x3 JT = transpose(J);
  JT = transpose(Q)*JT;
  return transpose(JT);
}

Vector3 computeEdgeTangentAtVertex(Mesh* m, MeshEntity* edge,
    MeshEntity* vert,
    const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(edge) == Mesh::EDGE);
  PCU_ALWAYS_ASSERT(m->getType(vert) == Mesh::VERTEX);

  // get the shared vertex
  MeshEntity* ev[2];

  m->getDownward(edge, 0, ev);

  PCU_ALWAYS_ASSERT(ev[0] == vert || ev[1] == vert);

  int index;
  if (ev[0] == vert)
    index = 0;
  else
    index = 1;

  Vector3 xi;
  if (index == 0)
    xi = Vector3(-1.0, 0.0, 0.0);
  else
    xi = Vector3( 1.0, 0.0, 0.0);

  MeshElement* edgeElem = createMeshElement(m, edge);
  Matrix3x3 J;
  getJacobian(edgeElem, xi, J);
  destroyMeshElement(edgeElem);

  J = getJacobianInMetric(J, Q);

  Vector3 t;
  if (index == 0)
    t = Vector3(J[0]);
  else
    t = Vector3(0.0, 0.0, 0.0) - Vector3(J[0]);
  return t.normalize();
}


Vector3 computeFaceNormalAtVertex(Mesh* m, MeshEntity* face,
    MeshEntity* vert,
    const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(face) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(vert) == Mesh::VERTEX);

  // get the index of the vertex in face
  MeshEntity* fv[3];
  m->getDownward(face, 0, fv);

  int index = findIn(fv, 3, vert);

  PCU_ALWAYS_ASSERT(index > -1 && index < 3);

  Vector3 xi;
  if (index == 0)
    xi = Vector3(0.0, 0.0, 0.0);
  else if (index == 1)
    xi = Vector3(1.0, 0.0, 0.0);
  else
    xi = Vector3(0.0, 1.0, 0.0);

  MeshElement* faceElem = createMeshElement(m, face);
  Matrix3x3 J;
  getJacobian(faceElem, xi, J);
  destroyMeshElement(faceElem);

  J = getJacobianInMetric(J, Q);

  return cross(J[0], J[1]).normalize();
}

static double computeEdgeEdgeCosAngleInTri(Mesh* m, MeshEntity* tri,
    MeshEntity* e1, MeshEntity* e2, const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tri) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(e1) == Mesh::EDGE);
  PCU_ALWAYS_ASSERT(m->getType(e2) == Mesh::EDGE);

  // get the shared vertex
  MeshEntity* e1v[2];
  MeshEntity* e2v[2];
  MeshEntity* sharedVert;

  m->getDownward(e1, 0, e1v);
  m->getDownward(e2, 0, e2v);

  if (e1v[0] == e2v[0] || e1v[0] == e2v[1])
    sharedVert = e1v[0];
  else
    sharedVert = e1v[1];

  PCU_ALWAYS_ASSERT(sharedVert);


  Vector3 t1 = computeEdgeTangentAtVertex(m, e1, sharedVert, Q);
  Vector3 t2 = computeEdgeTangentAtVertex(m, e2, sharedVert, Q);
  return t1*t2;
}

Vector3 computeFaceNormalAtEdgeInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* face, MeshEntity* edge, Matrix3x3 Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tet) == Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(face) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(edge) == Mesh::EDGE);

  MeshEntity* te[6];
  m->getDownward(tet, 1, te);
  int index = -1;
  for (int i = 0; i < 6; i++)
    if (te[i] == edge) {
      index = i;
      break;
    }
  PCU_ALWAYS_ASSERT(index > -1 && index < 6);
  Vector3 xi; // corresponding parametric coord in tet for the mid-edge point
  switch (index) {
    case 0:
      xi = Vector3(0.5, 0.0, 0.0);
      break;
    case 1:
      xi = Vector3(0.5, 0.5, 0.0);
      break;
    case 2:
      xi = Vector3(0.0, 0.5, 0.0);
      break;
    case 3:
      xi = Vector3(0.0, 0.0, 0.5);
      break;
    case 4:
      xi = Vector3(0.5, 0.0, 0.5);
      break;
    case 5:
      xi = Vector3(0.0, 0.5, 0.5);
      break;
    default:
      break;
  }
  MeshElement* tetElem = createMeshElement(m, tet);
  Matrix3x3 J;
  getJacobian(tetElem, xi, J);
  destroyMeshElement(tetElem);

  // get the two columns of Jacobian corresponding to face tangents
  MeshEntity* tf[4];
  m->getDownward(tet, 2, tf);
  index = -1;
  for (int i = 0; i < 4; i++)
    if (tf[i] == face) {
      index = i;
      break;
    }
  PCU_ALWAYS_ASSERT(index > -1 && index < 4);

  // transform Jacobian into Metric space
  J = getJacobianInMetric(J, Q);
  if (index == 0)
    return cross(J[0], J[1]).normalize();
  else if (index == 1)
    return cross(J[2], J[0]).normalize();
  else if (index == 2)
    return cross(J[2]-J[0], J[1]-J[0]).normalize();
  else
    return cross(J[1], J[2]).normalize();
}


static double computeFaceFaceCosAngleInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* face1, MeshEntity* face2, const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tet) == Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(face1) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(face2) == Mesh::TRIANGLE);

  if (face1 == face2) return 1.0;

  MeshEntity* f1e[3];
  MeshEntity* f2e[3];
  MeshEntity* sharedEdge = 0;

  m->getDownward(face1, 1, f1e);
  m->getDownward(face2, 1, f2e);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (f1e[i] == f2e[j]) {
      	sharedEdge = f1e[i];
      	break;
      }

  PCU_ALWAYS_ASSERT(sharedEdge);
  Vector3 f1n = computeFaceNormalAtEdgeInTet(m, tet, face1, sharedEdge, Q);
  Vector3 f2n = computeFaceNormalAtEdgeInTet(m, tet, face2, sharedEdge, Q);
  // turn f2n so you get the inside angle
  f2n = Vector3(0.0, 0.0, 0.0) - f2n;
  return f1n * f2n;
}

static double computeEdgeFaceCosAngleInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* edge, MeshEntity* face, const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tet) == Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(edge) == Mesh::EDGE);
  PCU_ALWAYS_ASSERT(m->getType(face) == Mesh::TRIANGLE);

  MeshEntity* fe[3];
  m->getDownward(face, 1, fe);
  if (findIn(fe, 3, edge) > -1)
    return 1.0;


  // find the share vertex
  MeshEntity* ev[2];
  MeshEntity* fv[3];
  MeshEntity* sharedVert = 0;

  m->getDownward(edge, 0, ev);
  m->getDownward(face, 0, fv);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      if (ev[i] == fv[j]) {
      	sharedVert = ev[i];
      	break;
      }

  PCU_ALWAYS_ASSERT(sharedVert);

  Vector3 t = computeEdgeTangentAtVertex(m, edge, sharedVert, Q);
  Vector3 n = computeFaceNormalAtVertex(m, face, sharedVert, Q);
  n = Vector3(0.0, 0.0, 0.0) - n;

  double tmp = t*n;
  tmp = 1. - tmp*tmp;
  return std::sqrt(tmp);
}


static double computeEdgeEdgeCosAngleInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* edge1, MeshEntity* edge2, const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tet) == Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(edge1) == Mesh::EDGE);
  PCU_ALWAYS_ASSERT(m->getType(edge2) == Mesh::EDGE);

  if (edge1 == edge2) return 1.0;

  // find the share vertex
  MeshEntity* e1v[2];
  MeshEntity* e2v[2];
  MeshEntity* sharedVert = 0;

  m->getDownward(edge1, 0, e1v);
  m->getDownward(edge2, 0, e2v);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      if (e1v[i] == e2v[j]) {
      	sharedVert = e1v[i];
      	break;
      }

  PCU_ALWAYS_ASSERT(sharedVert);

  Vector3 t1 = computeEdgeTangentAtVertex(m, edge1, sharedVert, Q);
  Vector3 t2 = computeEdgeTangentAtVertex(m, edge2, sharedVert, Q);

  return t1 * t2;
}


double computeCosAngleInTri(Mesh* m, MeshEntity* tri,
    MeshEntity* e1, MeshEntity* e2, const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tri) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(e1) == Mesh::EDGE);
  PCU_ALWAYS_ASSERT(m->getType(e2) == Mesh::EDGE);

  return computeEdgeEdgeCosAngleInTri(m, tri, e1, e2, Q);
}


double computeCosAngleInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* e1, MeshEntity* e2, const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT(m->getType(tet) == Mesh::TET);
  PCU_ALWAYS_ASSERT(m->getType(e1) == Mesh::EDGE ||
      m->getType(e1) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(e2) == Mesh::EDGE ||
      m->getType(e2) == Mesh::TRIANGLE);
  if (m->getType(e1) == Mesh::TRIANGLE && m->getType(e2) == Mesh::TRIANGLE)
    return computeFaceFaceCosAngleInTet(m, tet, e1, e2, Q);
  else if (m->getType(e1) == Mesh::TRIANGLE) // e1 is a face
    return computeEdgeFaceCosAngleInTet(m, tet, e2 /*edge*/, e1 /*face*/, Q);
  else if (m->getType(e2) == Mesh::TRIANGLE) // e2 is a face
    return computeEdgeFaceCosAngleInTet(m, tet, e1 /*edge*/, e2 /*face*/, Q);
  else // both e1 and e2 are edges
    return computeEdgeEdgeCosAngleInTet(m, tet, e1, e2, Q);
}

static double getEdgeLength(Mesh* m, MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(m->getType(e) == Mesh::EDGE);
  MeshEntity* vs[2];
  m->getDownward(e, 0, vs);
  Vector3 p0, p1;
  m->getPoint(vs[0], 0, p0);
  m->getPoint(vs[1], 0, p1);
  return (p1 - p0).getLength();
}

double computeShortestHeightInTet(Mesh* m, MeshEntity* tet,
    const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT_VERBOSE(m->getType(tet) == Mesh::TET,
      "Expecting a tet. Aborting! ");

  double minHeight = 1.0e12;

  MeshEntity* faces[4];
  MeshEntity* edges[6];

  m->getDownward(tet, 2, faces);
  m->getDownward(tet, 1, edges);

  // iterate over faces
  for (int i = 0; i < 4; i++) {
    MeshEntity* currentFace = faces[i];
    MeshEntity* outOfFaceEdge = 0;
    MeshEntity* es[3];
    m->getDownward(currentFace, 1, es);
    // get the out of face edge
    for (int j = 0; j < 6; j++) {
      int index = findIn(es, 3, edges[j]);
      if (index == -1) {
      	outOfFaceEdge = edges[j];
      	break;
      }
    }
    PCU_ALWAYS_ASSERT(outOfFaceEdge);
    // compute the cos angle
    double edgeFaceCosAngle =
      computeEdgeFaceCosAngleInTet(m, tet, outOfFaceEdge, currentFace, Q);
    // compute the sin angle
    double edgeFaceSinAngle = std::sqrt(
    	1 - edgeFaceCosAngle*edgeFaceCosAngle);
    // height is sinAngle*edgeLen
    double currentHeight = edgeFaceSinAngle * getEdgeLength(m, outOfFaceEdge);

    if (currentHeight < minHeight)
      minHeight = currentHeight;
  }
  return minHeight;
}

double computeShortestHeightInTri(Mesh* m, MeshEntity* tri,
    const Matrix3x3& Q)
{
  PCU_ALWAYS_ASSERT_VERBOSE(m->getType(tri) == Mesh::TRIANGLE,
      "Expecting a tri. Aborting! ");

  double minHeight = 1.0e12;

  MeshEntity* edges[3];

  m->getDownward(tri, 1, edges);

  // iterate over edges
  for (int i = 0; i < 3; i++) {
    MeshEntity* currentEdge = edges[i];
    MeshEntity*	   nextEdge = edges[(i+1)%3];
    // compute the cos angle
    double edgeEdgeCosAngle =
      computeEdgeEdgeCosAngleInTri(m, tri, currentEdge, nextEdge, Q);
    // compute the sin angle
    double edgeEdgeSinAngle = std::sqrt(
    	1 - edgeEdgeCosAngle*edgeEdgeCosAngle);
    // height is sinAngle*edgeLen
    double currentHeight = edgeEdgeSinAngle * getEdgeLength(m, nextEdge);

    if (currentHeight < minHeight)
      minHeight = currentHeight;
  }
  return minHeight;
}

}
