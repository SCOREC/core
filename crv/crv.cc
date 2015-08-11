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
  static int table[13] = {1,1,2,6,24,120,720,5040,40320,362880,
      3628800,39916800,479001600};
  return table[i];
}

int binomial(int n, int i)
{

  i = std::min(n-i,i);

  if(i == 0)
    return 1;
  if(i == 1)
    return n;

  static int const bn4[1] = {6};
  static int const bn5[1] = {10};
  static int const bn6[2] = {15,20};
  static int const bn7[2] = {21,35};
  static int const bn8[3] = {28,56,70};
  static int const bn9[3] = {36,84,126};
  static int const bn10[4] = {45,120,210,252};
  static int const bn11[4] = {55,165,330,462};
  static int const bn12[5] = {66,220,495,792};
  static int const bn13[5] = {78,286,715,1287,1716};
  static int const* const bnTable[10] = {bn4,bn5,bn6,bn7,bn8,
      bn9,bn10,bn11,bn12,bn13};

  return bnTable[n-4][i-2];
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
  if (m->getModelType(g) == 3)
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
