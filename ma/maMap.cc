/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#include "maMap.h"
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
  static GetMapFunction table[TYPES] =
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
  if (type == EDGE)
    return std::min(xi[0] + 1, 1 - xi[0]);
/* returns the least barycentric coordinate */
  if (type == TRI)
    return std::min(xi[0],
           std::min(xi[1],
                    1 - xi[0] - xi[1]));
  if (type == TET)
    return std::min(xi[0],
           std::min(xi[1],
           std::min(xi[2],
                    1 - xi[0] - xi[1] - xi[2])));
  return 0;
}

}
