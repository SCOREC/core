/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#include "maMap.h"
#include <pcu_util.h>
#include <algorithm>

namespace ma {

void getVertPoints(apf::Mesh* m, Entity* e, Vector* p)
{
  apf::Downward v;
  int n = m->getDownward(e,0,v);
  for (int i = 0; i < n; ++i)
    m->getPoint(v[i],0,p[i]);
}

/* dealing with our choice of using [-1,1]
   for edge coordinates... */
Affine getEdgeMap(apf::Mesh* m, Entity* e)
{
  Vector p[2];
  getVertPoints(m,e,p);
  Vector mp = (p[1] + p[0]) / 2;
  Affine a;
  a.A = apf::getFrame(p[1] - mp);
  a.A = apf::transpose(a.A);
  a.b = mp;
  return a;
}

Affine getTriMap(apf::Mesh* m, Entity* e)
{
  Vector p[3];
  getVertPoints(m,e,p);
  Affine a;
  a.A[0] = p[1] - p[0];
  a.A[1] = p[2] - p[0];
  a.A[2] = apf::cross(a.A[0],a.A[1]);
  a.A = apf::transpose(a.A);
  a.b = p[0];
  return a;
}

Affine getTetMap(apf::Mesh* m, Entity* e)
{
  Vector p[4];
  getVertPoints(m,e,p);
  Affine a;
  a.A[0] = p[1] - p[0];
  a.A[1] = p[2] - p[0];
  a.A[2] = p[3] - p[0];
  a.A = apf::transpose(a.A);
  a.b = p[0];
  return a;
}

Affine getMap(apf::Mesh* m, Entity* e)
{
  typedef Affine (*GetMapFunction)(apf::Mesh*,Entity*);
  static GetMapFunction table[apf::Mesh::TYPES] =
  {0,
   getEdgeMap,
   getTriMap,
   0,
   getTetMap,
   0,
   0,
   0};
  return table[m->getType(e)](m,e);
}

/* negative values are outside the parent domain,
   positive values are inside.
   increasing distance from the centroid decreases this value,
   but using a simple vector distance from the centroid does
   not have the required property that all inside values are greater
   than all outside values */
double getInsideness(apf::Mesh* m, Entity* e, Vector const& xi)
{
  int type = m->getType(e);
  if (type == apf::Mesh::EDGE)
    return std::min(xi[0] + 1, 1 - xi[0]);
/* returns the least barycentric coordinate */
  if (type == apf::Mesh::TRIANGLE)
    return std::min(xi[0],
           std::min(xi[1],
                    1 - xi[0] - xi[1]));
  if (type == apf::Mesh::TET)
    return std::min(xi[0],
           std::min(xi[1],
           std::min(xi[2],
                    1 - xi[0] - xi[1] - xi[2])));
  return 0;
}

Vector curvedElemInvMap(
    apf::Mesh* m,
    Entity* e,
    const Vector& p,
    const double tol,
    const int maxIter)
{
  int iter = 0;
  // put the initial guess at the center of the elements
  Vector xi0;
  int type = m->getType(e);
  if (type == apf::Mesh::VERTEX)
    xi0 = Vector(0., 0., 0.);
  else if (type == apf::Mesh::EDGE)
    xi0 = Vector(0., 0., 0.);
  else if (type == apf::Mesh::TRIANGLE)
    xi0 = Vector(1./3., 1./3., 1./3.);
  else if (type == apf::Mesh::TET)
    xi0 = Vector(1./4., 1./4., 1./4.);
  else
    PCU_ALWAYS_ASSERT_VERBOSE(0, "unsupported type!");


  apf::MeshElement* me = apf::createMeshElement(m, e);

  Matrix Jinv;
  Vector x;
  Vector xi_new = xi0;
  Vector xi_old = xi0;
  double err = 1.e16;

  // initial iteration
  apf::getJacobianInv(me, xi_old, Jinv);
  Jinv = apf::transpose(Jinv);
  apf::mapLocalToGlobal(me, xi_old, x);
  xi_new = xi_old - Jinv*(x-p);
  err = (xi_new - xi_old)*(xi_new - xi_old);

  while (err > tol) {
    iter++;
    if (iter > maxIter) break;
    xi_old = xi_new;
    apf::getJacobianInv(me, xi_old, Jinv);
    Jinv = apf::transpose(Jinv);
    apf::mapLocalToGlobal(me, xi_old, x);
    xi_new = xi_old - Jinv*(x-p);
    err = (xi_new - xi_old)*(xi_new - xi_old);
  }

  apf::destroyMeshElement(me);
  return xi_new;
}

}
