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

int getNumInternalControlPoints(int type, int order)
{
  switch (type) {
    case apf::Mesh::EDGE:
      return order-1;
      break;
    case apf::Mesh::TRIANGLE:
      return (order-1)*(order-2)/2;
      break;
    case apf::Mesh::TET:
      return (order-1)*(order-2)*(order-3)/6;
      break;
    default:
      break;
  }
  return 0;
}

int getNumControlPoints(int type, int order)
{
  switch (type) {
    case apf::Mesh::EDGE:
      return order+1;
      break;
    case apf::Mesh::TRIANGLE:
      return (order+1)*(order+2)/2;
      break;
    case apf::Mesh::TET:
      return (order+1)*(order+2)*(order+3)/6;
      break;
    default:
      break;
  }
  return 0;
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
