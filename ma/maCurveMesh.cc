/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/
#include "maCurveMesh.h"
#include "maSnap.h"
#include <apfField.h>
#include <apfShape.h>
#include <apfMesh.h>
#include <gmi.h>

#include <string.h>
#include <stdio.h> //TODO: dont need this normally
#include <sstream>

namespace ma {

double interpolationError(Mesh* m, Entity* e, int n,
    Vector &samplept, Vector &maxpt){
  Model* g = m->toModel(e);
  if (m->getModelType(g) == m->getDimension()) return 0.;
  int d = apf::getDimension(m,e);
  int nj = (d == 2) ? n : 1;
  Vector pt,pa(0.,0.,0.),cpt,cpa;
  double max = -1.0;
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  for (int j = 0; j < nj; ++j){
    pa[1] = 1.*j/(nj-1.);
    for (int i = 0; i < n-j; ++i){
      if(d == 1) pa[0] = 2.*i/(n-1)-1.0;
      else pa[0] = 1.*i/(nj-1.);
      apf::getVector(elem,pa,pt);
      m->getClosestPoint(g,pt,cpt,cpa);
      if((cpt-pt).getLength() > max){
        max = (cpt-pt).getLength();
        maxpt = cpt;
        samplept = pt;
      }
    }
  }

  apf::destroyElement(elem);
  return max;
}

double interpolationErrorAtNodeXi(Mesh* m, Entity* e, 
    Vector &samplept, Vector &maxpt)
{
  Model* g = m->toModel(e);
  if (m->getModelType(g) == m->getDimension())
    return 0.;
  int d = apf::getDimension(m,e);
  Vector pt,pa(0.,0.,0.),cpt,cpa;
  double max = -1.0;
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  apf::FieldShape * fs = m->getCoordinateField()->getShape();
  for (int i = 0; i < fs->countNodesOn(d); ++i){
    fs->getNodeXi(d,i,pa);
    apf::getVector(elem,pa,pt);
    m->getClosestPoint(g,pt,cpt,cpa);
    if((cpt-pt).getLength() > max){
      max = (cpt-pt).getLength();
      maxpt = cpt;
      samplept = pt;
    }
  }

  apf::destroyElement(elem);
  return max;
}

static void snapToInterpolate(Mesh* m, int d){

  apf::FieldShape * fs = m->getCoordinateField()->getShape();
  int t = (d == 1) ? Mesh::EDGE : Mesh::TRIANGLE;
  int non = fs->countNodesOn(t);
  Entity* e;
  Vector p, xi, pt;
  Iterator* it = m->begin(d);
  while ((e = m->iterate(it))) {
    Model* g = m->toModel(e);
    if(m->getModelType(g) == m->getDimension()) continue;
    for(int i = 0; i < non; ++i){
      fs->getNodeXi(t,i,xi);
      if(t == Mesh::EDGE)
        transferParametricOnEdgeSplit(m,e,0.5*(xi[0]+1.),p);
      else
        transferParametricOnTriSplit(m,e,xi,p);
      m->snapToModel(g,p,pt);
      m->setPoint(e,i,pt);
    }
  }
  m->end(it);
}

void convertInterpToBezier(Mesh* m, Entity* e, int n, int ne,
    apf::NewArray<double>& c){

  apf::NewArray<Vector> l(n), b(ne);
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  apf::getVectorNodes(elem,l);

  for(int i = 0; i < ne; ++i)
    b[i] = Vector(0,0,0);

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      b[i] += l[j]*c[i*n+j];

  for(int i = 0; i < ne; ++i)
    m->setPoint(e,i,b[i]);

  apf::destroyElement(elem);
}

void curveMeshToBezier(Mesh* m, int order){
  assert(order < 7);
  assert(order > 0);

  int md = m->getDimension();

  apf::changeMeshShape(m, apf::getBezier(md,order),true);

  apf::FieldShape * fs = m->getCoordinateField()->getShape();
  printf("Changing Mesh Shape to %s\n", fs->getName());
  if(order == 1) return;

  Entity* e;

  // interpolate points in each dimension
  for(int d = 1; d < md; ++d)
    snapToInterpolate(m,d);

  // go downward, and convert interpolating to control points
  for(int d = md-1; d >= 1; --d){
    int n = (d == 2)? (order+1)*(order+2)/2 : order+1;
    int ne = fs->countNodesOn(d);

    apf::NewArray<double> c;
    apf::getTransformationCoefficients(order,md,d,c);

    Iterator* it = m->begin(d);
    while ((e = m->iterate(it))) {
      Model* g = m->toModel(e);
      if(m->getModelType(g) == m->getDimension()) continue;

      convertInterpToBezier(m,e,n,ne,c);

    }
    m->end(it);
  }

  printf("Done Changing Mesh Shape\n");
  m->acceptChanges();
  m->verify();
}


void writePointSet(Mesh* m, int d, int n, const char* prefix){
  int nj = (d == 2) ? n : 1;
  apf::DynamicArray<apf::Vector3> pts(0);
  Iterator* it = m->begin(d);
  Entity* e;
  Vector pa(0.,0.,0.),pt(0.,0.,0.);
  while ((e = m->iterate(it))) {

    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);

    for (int j = 0; j < nj+1; ++j){
      pa[1] = 1.*j/nj;
      for (int i = 0; i < n+1-j; ++i){
        if(d == 1) pa[0] = 2.*i/n-1.;
        else pa[0] = 1.*i/nj;
        apf::getVector(elem,pa,pt);
        pts.append(pt);
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);
  std::stringstream ss;
  ss << prefix << "_" << d << "_"
      << m->getCoordinateField()->getShape()->getOrder();
  apf::writeCSVPointSet(ss.str().c_str(),pts);
}

} // namespace ma

