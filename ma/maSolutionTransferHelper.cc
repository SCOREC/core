
/*******************************************************************************

opyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSolutionTransferHelper.h"
#include "maAffine.h"
#include "maMap.h"
#include <apfShape.h>
#include <apfNumbering.h>
#include <float.h>

namespace ma {

int getMinimumDimension(apf::FieldShape* s)
{
  int transferDimension = 4;
  for (int d=1; d <= 3; ++d)
    if (s->hasNodesIn(d))
    {
      transferDimension = d;
      break;
    }
  return transferDimension;
}

int getBestElement(
    apf::Mesh* mesh,
    int n,
    apf::Element** elems,
    Affine* elemInvMaps,
    Vector const& point,
    Vector& bestXi)
{
  double bestValue = -DBL_MAX;
  int bestI = 0;
  for (int i = 0; i < n; ++i)
  {
    Vector xi;
    if (mesh->getShape()->getOrder() == 1)
      xi = elemInvMaps[i] * point;
    else
      xi = curvedElemInvMap(mesh, apf::getMeshEntity(elems[i]), point);
    double value = getInsideness(mesh,apf::getMeshEntity(elems[i]),xi);
    if (value > bestValue)
    {
      bestValue = value;
      bestI = i;
      bestXi = xi;
    }
  }
  return bestI;
}

void transferToNode(
    apf::Field* field,
    double *value,
    int n,
    apf::Element** elems,
    Affine* elemInvMaps,
    apf::Node const& node)
{
  apf::Mesh* mesh = apf::getMesh(field);
  // first get the physical coordinate of the node
  Vector xi;
  Vector point;
  apf::FieldShape* shape = apf::getShape(field);
  shape->getNodeXi(mesh->getType(node.entity),node.node,xi);
  if (mesh->getShape()->getOrder() == 1) { // if linear mesh use the affine
    Affine childMap = getMap(mesh,node.entity);
    point = childMap * xi;
  }
  else { // else inquire the physical coordinated of local coordinate xi
    apf::MeshElement* me = apf::createMeshElement(mesh,node.entity);
    apf::mapLocalToGlobal(me, xi, point);
    apf::destroyMeshElement(me);
  }
  Vector elemXi;
  int i = getBestElement(mesh,n,elems,elemInvMaps,point,elemXi);
  apf::getComponents(elems[i],elemXi,value);
  apf::setComponents(field,node.entity,node.node,value);
}

void transfer(
    apf::Field* field,
    double *value,
    int n, // size of the cavity
    Entity** cavity,
    EntityArray& newEntities)
{
  apf::Mesh* mesh = apf::getMesh(field);
  apf::FieldShape* shape = apf::getShape(field);
  int minDim = getMinimumDimension(shape);
  if (getDimension(mesh, cavity[0]) < minDim)
    return;
  apf::NewArray<apf::Element*> elems(n);
  for (int i = 0; i < n; ++i)
    elems[i] = apf::createElement(field,cavity[i]);
  apf::NewArray<Affine> elemInvMaps(n);
  for (int i = 0; i < n; ++i)
    elemInvMaps[i] = invert(getMap(mesh,cavity[i]));
  for (size_t i = 0; i < newEntities.getSize(); ++i)
  {
    int type = mesh->getType(newEntities[i]);
    if (type == apf::Mesh::VERTEX)
      continue; //vertices will have been handled specially beforehand
    int nnodes = shape->countNodesOn(type);
    for (int j = 0; j < nnodes; ++j)
    {
      apf::Node node(newEntities[i],j);
      transferToNode(field, value,
      	  n,&(elems[0]),&(elemInvMaps[0]),node);
    }
  }
  for (int i = 0; i < n; ++i)
    apf::destroyElement(elems[i]);
}
}
