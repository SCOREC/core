/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"

namespace crv {

int factorial(int i)
{
  static int table[11] = {1,1,2,6,24,120,720,5040,40320,362880,3628800};
  return table[i];
}

int binomial(int n, int i)
{
  static int const table[28] =
  {1,1,1,1,1,1,1,1,2,3,4,5,6,1,3,6,10,15,1,4,10,20,1,5,15,1,6,1};
  return table[i*7 - (i-1)*i/2 + n-i];
}

int trinomial(int n, int i, int j)
{
  return binomial(n,i)*binomial(n-i,j);
}

int quadnomial(int n, int i, int j, int k)
{
  return binomial(n,i)*binomial(n-i,j)*binomial(n-i-j,k);
}

void elevateBezierCurve(apf::Mesh2* m, apf::MeshEntity* edge, int n, int r)
{
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),edge);

  apf::Vector3 pt;
  apf::NewArray<apf::Vector3> p;
  apf::getVectorNodes(elem,p);

  assert(m->getType(edge) == apf::Mesh::EDGE);

  // reorder p into geometric ordering
  apf::NewArray<int> map(n+1);
  map[0] = 0; map[n] = 1;

  for(int i = 1; i < n; ++i)
    map[i] = i+1;
  for(int i = 1; i < n+r; ++i){
    pt.zero();
    for(int j = std::max(0,i-r); j <= std::min(i,n); ++j)
      pt += p[map[j]]*binomial(n,j)*
      binomial(r,i-j)/binomial(n+r,i);
    m->setPoint(edge,i-1,pt);
  }

  apf::destroyElement(elem);
}

double interpolationError(apf::Mesh* m, apf::MeshEntity* e, int n){
  apf::ModelEntity* g = m->toModel(e);
  if (m->getModelType(g) == m->getDimension())
    return 0.;
  int d = apf::getDimension(m,e);
  int nj = (d == 2) ? n : 1;
  apf::Vector3 pt,pa(0.,0.,0.),cpt,cpa;
  double max = 0.0;
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  for (int j = 0; j <= nj; ++j){
    pa[1] = 1.*j/nj;
    for (int i = 0; i <= n-j; ++i){
      if(d == 1) pa[0] = 2.*i/n-1.;
      else pa[0] = 1.*i/n;
      apf::getVector(elem,pa,pt);
      m->getClosestPoint(g,pt,cpt,cpa);
      max = std::max((cpt-pt).getLength(),max);
    }
  }
  apf::destroyElement(elem);
  return max;
}

void fail(const char* why)
{
  fprintf(stderr,"CRV FAILED: %s\n",why);
  abort();
}

}
