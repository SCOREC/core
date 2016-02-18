/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvTables.h"
#include <cstdlib>
#include <iostream>

namespace crv {

int getNumInternalControlPoints(int type, int order)
{
  assert(order > 0);
  switch (type) {
    case apf::Mesh::VERTEX:
      return 1;
      break;
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
  fail("invalid type/order combination\n");
  return 0;
}

int getNumControlPoints(int type, int order)
{
  assert(order > 0);
  switch (type) {
    case apf::Mesh::VERTEX:
      return 1;
      break;
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
  fail("invalid type/order combination\n");
  return 0;
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

void getTransformationMatrix(apf::Mesh* m, apf::MeshEntity* e,
    mth::Matrix<double>& A, const apf::Vector3 *range)
{

  int type = m->getType(e);
  apf::FieldShape* fs = m->getShape();
  apf::EntityShape* es = fs->getEntityShape(type);
  int n = es->countNodes();
  int typeDim = apf::Mesh::typeDimension[type];

  apf::Vector3 xi, exi;
  int evi = 0;
  apf::NewArray<double> values;
  apf::NewArray<double> shape_vals;

  A.zero();

  int row = 0;
  // loop over lower entities to get the xi at each point
  for(int d = 0; d <= typeDim; ++d){
    int nDown = apf::Mesh::adjacentCount[type][d];
    int bt = apf::Mesh::simplexTypes[d];
    apf::EntityShape* shape = apf::getLagrange(1)->getEntityShape(bt);
    int non = fs->countNodesOn(bt);
    int nvtx =  apf::Mesh::adjacentCount[bt][0];
    for(int j = 0; j < nDown; ++j){
      for(int x = 0; x < non; ++x){
        fs->getNodeXi(bt,x,xi);
        shape->getValues(0, 0, xi, shape_vals);
        // get the xi of the lower entity wrt to the main entity
        if(d < typeDim){
          exi.zero();
          evi = j;
          for (int i = 0; i < nvtx; ++i) {
            if(bt == apf::Mesh::EDGE && type == apf::Mesh::TRIANGLE)
              evi = apf::tri_edge_verts[j][i];
            if(bt == apf::Mesh::EDGE && type == apf::Mesh::TET)
              evi = apf::tet_edge_verts[j][i];
            if(bt == apf::Mesh::TRIANGLE && type == apf::Mesh::TET)
              evi = apf::tet_tri_verts[j][i];
            exi += elem_vert_xi[type][evi] * shape_vals[i];
          }
        } else {
          exi = xi;
        }
        // scale to the required range
        // slight change here because edges run [-1,1]
        if(typeDim == 1){
          exi[0] = 0.5*(exi[0]+1);
          exi[0] = range[1][0]*exi[0]+range[0][0]*(1.-exi[0]);
        }
        if(typeDim == 2){
          double p = 1.-exi[0]-exi[1];
          exi = range[0]*p + range[1]*exi[0]+range[2]*exi[1];
        }
        if(typeDim == 3){
          double p = 1.-exi[0]-exi[1]-exi[2];
          exi = range[0]*p + range[1]*exi[0]+range[2]*exi[1]+range[3]*exi[2];
        }
        es->getValues(m,e,exi,values);
        for(int i = 0; i < n; ++i){
          A(row,i) = values[i];
        }
        ++row;
      }
    }
  }
}

int countNumberInvalidElements(apf::Mesh2* m)
{
  int n = 0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(m->getDimension());
  if (m->getShape()->getOrder() == 1){
    while ((e = m->iterate(it))) {
      n += (apf::measure(m,e) < 1e-10);
    }
  } else {
    while ((e = m->iterate(it))) {
      n += (checkValidity(m,e,4) > 1);
    }
  }
  m->end(it);
  return n;
}

void fail(const char* why)
{
  fprintf(stderr,"CRV FAILED: %s\n",why);
  abort();
}

}
