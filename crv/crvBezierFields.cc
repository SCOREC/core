/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"
#include "crvShape.h"
#include "crvTables.h"
#include "crvQuality.h"
#include <apfField.h>
#include <lionPrint.h>
#include <cstdlib>
#include <mth_def.h>
#include <iostream>
#include <pcu_util.h>

namespace crv {

static void convertVectorFieldInterpolationPoints(int n, int ne,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<double>& c,
    apf::NewArray<apf::Vector3>& newNodes){

  for(int i = 0; i < ne; ++i)
    newNodes[i].zero();

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      newNodes[i] += nodes[j]*c[i*n+j];
}

static void convertScalarFieldInterpolationPoints(int n, int ne,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& c,
    apf::NewArray<double>& newNodes)
{

  for(int i = 0; i < ne; ++i)
    newNodes[i] = 0.;

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      newNodes[i] += nodes[j]*c[i*n+j];
}

void convertInterpolationFieldPoints(apf::MeshEntity* e,
    apf::Field* f,
    int n, int ne, apf::NewArray<double>& c)
{

  apf::NewArray<apf::Vector3> l, b(ne);
  apf::NewArray<double> ls, bs(ne);
  apf::MeshElement* me =
      apf::createMeshElement(f->getMesh(),e);
  apf::Element* elem =
      apf::createElement(f,me);

  if (apf::getValueType(f) == apf::VECTOR) {
    apf::getVectorNodes(elem,l);
    convertVectorFieldInterpolationPoints(n, ne, l, c, b);
    for(int i = 0; i < ne; ++i)
      apf::setVector(f, e, i, b[i]);
  }
  else if (apf::getValueType(f) == apf::SCALAR) {
    apf::getScalarNodes(elem, ls);
    convertScalarFieldInterpolationPoints(n, ne, ls, c, bs);
    for(int i = 0; i < ne; ++i)
      apf::setScalar(f, e, i, bs[i]);
  }
  else
    printf("Field type not implemented\n");
  apf::destroyElement(elem);
  apf::destroyMeshElement(me);
}

void convertInterpolatingFieldToBezier(apf::Mesh2* m_mesh, apf::Field* f)
{
  // TODO: to be completed
  // Zero out the blending orders (since they are used internally
  // in Bezier shapes (triangles and tets) but save the old ones to
  // revert back to at the end of this (so the rest of Bezier related
  // procedures work as intendedi)
  //
  //setBlendingOrder(apf::Mesh::TYPES, 0);

  apf::FieldShape * fs = apf::getShape(f);
  int order = fs->getOrder();

  //apf::Field* fnew = createField(
  //    m_mesh, "copy_field", apf::getValueType(f), crv::getBezier(order));
  //transferFields(m_mesh, f, fnew);

  int md = m_mesh->getDimension();
  //
  int blendingOrder = getBlendingOrder(apf::Mesh::simplexTypes[md]);
  //blendingOrder = 0;
  // go downward, and convert interpolating to control points
  int startDim = md - (blendingOrder > 0);

  for(int d = startDim; d >= 1; --d){
  //for(int d = md; d >= 1; d--){
    if(!fs->hasNodesIn(d)) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<double> c;
    getBezierTransformationCoefficients(order,
        apf::Mesh::simplexTypes[d],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      if(m_mesh->isOwned(e))
        convertInterpolationFieldPoints(e,f,n,ne,c);
    }
    m_mesh->end(it);
  }

  for( int d = 2; d <= md; ++d){
    if(!fs->hasNodesIn(d) ||
        getBlendingOrder(apf::Mesh::simplexTypes[d])) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<double> c;
    getInternalBezierTransformationCoefficients(m_mesh,order,1,
        apf::Mesh::simplexTypes[d],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      if(!isBoundaryEntity(m_mesh,e) && m_mesh->isOwned(e))
        convertInterpolationFieldPoints(e,f,n-ne,ne,c);
    }
    m_mesh->end(it);
  }

  apf::synchronize(f);

}

}
