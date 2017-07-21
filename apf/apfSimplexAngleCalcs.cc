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

static Vector3 computeEdgeTangentAtVertex(Mesh* m, MeshEntity* edge,
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

static Vector3 computeFaceNormalAtEdge(Mesh* m, /*MeshEntity* tet,*/
    MeshEntity* face, MeshEntity* edge, Matrix3x3 Q)
{
  PCU_ALWAYS_ASSERT(m->getType(face) == Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(edge) == Mesh::EDGE);

  MeshEntity* fe[3];
  m->getDownward(face, 1, fe);
  int index = -1;
  for (int i = 0; i < 3; i++)
    if (fe[i] == edge) {
      index = i;
      break;
    }
  PCU_ALWAYS_ASSERT(index > -1 && index < 3);
  Vector3 xi; // corresponding parametric coord in tri of the mid-edge point
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
    default:
      break;
  }
  MeshElement* faceElem = createMeshElement(m, face);
  Matrix3x3 J;
  getJacobian(faceElem, xi, J);
  destroyMeshElement(faceElem);

  // transform Jacobian into Metric space
  J = getJacobianInMetric(J, Q);
  return cross(J[0], J[1]).normalize();
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
  Vector3 f1n = computeFaceNormalAtEdge(m, /*tet,*/ face1, sharedEdge, Q);
  Vector3 f2n = computeFaceNormalAtEdge(m, /*tet,*/ face2, sharedEdge, Q);
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

}
