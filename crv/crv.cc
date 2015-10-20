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

void getTransformationMatrix(apf::Mesh* m, apf::MeshEntity* e,
    mth::Matrix<double>& A)
{

  apf::Vector3 const edge_vert_xi[2] = {
      apf::Vector3(-1,0,0),
      apf::Vector3(1,0,0),
  };
  apf::Vector3 const tri_vert_xi[3] = {
      apf::Vector3(0,0,0),
      apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),
  };
  apf::Vector3 const tet_vert_xi[4] = {
      apf::Vector3(0,0,0),
      apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),
      apf::Vector3(0,0,1),
  };
  apf::Vector3 const* const elem_vert_xi[apf::Mesh::TYPES] = {
      0, /* vertex */
      edge_vert_xi,
      tri_vert_xi,
      0, /* quad */
      tet_vert_xi,
      0, /* hex */
      0, /* prism */
      0  /* pyramid */
  };

  int type = m->getType(e);
  apf::FieldShape* fs = m->getShape();
  apf::EntityShape* es = fs->getEntityShape(type);
  int n = es->countNodes();
  int typeDim = apf::Mesh::typeDimension[type];

  apf::Vector3 xi, exi;
  int evi = 0;
  apf::NewArray<double> values;

  A.zero();

  int boundaryTypes[4] = {apf::Mesh::VERTEX,apf::Mesh::EDGE,
      apf::Mesh::TRIANGLE,apf::Mesh::TET};

  int row = 0;
  for(int d = 0; d <= typeDim; ++d){
    int nDown = apf::Mesh::adjacentCount[type][d];
    for(int j = 0; j < nDown; ++j){
      int bt = boundaryTypes[d];
      apf::EntityShape* shape = apf::getLagrange(1)->getEntityShape(bt);

      for(int x = 0; x < fs->countNodesOn(bt); ++x){
        fs->getNodeXi(bt,x,xi);
        apf::NewArray<double> shape_vals;
        shape->getValues(0, 0, xi, shape_vals);

        if(d < typeDim){
          exi.zero();
          evi = j;
          for (int i = 0; i < apf::Mesh::adjacentCount[bt][0]; ++i) {
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
        es->getValues(m,e,exi,values);
        for(int i = 0; i < n; ++i){
          A(row,i) = values[i];
        }
        ++row;
      }
    }
  }
}

void fail(const char* why)
{
  fprintf(stderr,"CRV FAILED: %s\n",why);
  abort();
}

}
