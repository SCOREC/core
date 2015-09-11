/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include <cstdlib>

namespace crv {

/*
 * 18 is the maximum in the table, given that for n > 18,
 * quadnomial(n,i,j,k) can exceed MAX_INT and long's would be needed
 * This is also an upper bound on the order of tets, and implies a max order
 * of 7 to guarantee the full bezier jacobian determinant can work
 */
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
  static int const bn14[6] = {91,364,1001,2002,3003,3432};
  static int const bn15[6] = {105,455,1365,3003,5005,6435};
  static int const bn16[7] = {120,560,1820,4368,8008,11440,12870};
  static int const bn17[7] = {136,680,2380,6188,12376,19448,24310};
  static int const bn18[8] = {153,816,3060,8568,18564,31824,43758,48620};

  static int const* const bnTable[19] = {0,0,0,0,bn4,bn5,bn6,bn7,bn8,
      bn9,bn10,bn11,bn12,bn13,bn14,bn15,bn16,bn17,bn18};

  return bnTable[n][i-2];

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
  apf::NewArray<apf::Vector3> nodes, elevatedNodes(n+r+1);
  apf::getVectorNodes(elem,nodes);
  elevateBezierEdge(n,r,nodes,elevatedNodes);

  for(int i = 1; i < n+r; ++i)
    m->setPoint(edge,i-1,elevatedNodes[i]);

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
