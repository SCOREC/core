/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSolutionTransfer.h"
#include "maAffine.h"
#include "maMap.h"
#include "apfField.h"
#include <iostream>
namespace ma {

SolutionTransfer::~SolutionTransfer()
{
}

void SolutionTransfer::onVertex(
        apf::MeshElement*,
        Vector const&, 
        Entity*)
{
}

void SolutionTransfer::onRefine(
    Entity*,
    EntityArray&)
{
}

void SolutionTransfer::onCavity(
    EntityArray&,
    EntityArray&)
{
}

static int getMinimumDimension(apf::FieldShape* s)
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

int SolutionTransfer::getTransferDimension()
{
  int transferDimension = 4;
  for (int d=1; d <= 3; ++d)
    if (hasNodesOn(d))
    {
      transferDimension = d;
      break;
    }
  return transferDimension;
}
FieldTransfer::FieldTransfer(apf::Field* f)
{
  field = f;
  mesh = apf::getMesh(f);
  shape = apf::getShape(f);
  value.allocate(apf::countComponents(f));
}

bool FieldTransfer::hasNodesOn(int dimension)
{
  return shape->hasNodesIn(dimension);
}
void LinearTransfer::onVertex(
    apf::MeshElement* parent,
    Vector const& xi,
    Entity* vert)
{
  apf::Element* e = apf::createElement(field,parent);
  apf::getComponents(e,xi,&(value[0]));
  apf::setComponents(field,vert,0,&(value[0]));
  apf::destroyElement(e);
}

CavityTransfer::CavityTransfer(apf::Field* f):
  FieldTransfer(f)
{
  minDim = getMinimumDimension(getShape(f));
}
void CavityTransfer::transferToNodeIn(
    apf::Element* elem,
    apf::Node const& node,
    Vector const& elemXi)
{
  apf::getComponents(elem,elemXi,&(value[0]));
  apf::setComponents(field,node.entity,node.node,&(value[0]));
}
int CavityTransfer::getBestElement(
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
    Vector xi = elemInvMaps[i] * point;
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
void CavityTransfer::transferToNode(
    int n,
    apf::Element** elems,
    Affine* elemInvMaps,
    apf::Node const& node)
{
  Vector xi;
  shape->getNodeXi(mesh->getType(node.entity),node.node,xi);
  Affine childMap = getMap(mesh,node.entity);
  Vector point = childMap * xi;
  Vector elemXi;
  int i = getBestElement(n,elems,elemInvMaps,point,elemXi);
  transferToNodeIn(elems[i],node,elemXi);
}
void CavityTransfer::transfer(
    int n,
    Entity** cavity,
    EntityArray& newEntities)
{
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
      transferToNode(n,&(elems[0]),&(elemInvMaps[0]),node);
    }
  }
  for (int i = 0; i < n; ++i)
    apf::destroyElement(elems[i]);
}
void CavityTransfer::onRefine(
    Entity* parent,
    EntityArray& newEntities)
{
  apf::FieldShape *fs = apf::getShape(field);
  std::string name = fs->getName();
  if (name != std::string("Constant_3"))
    transfer(1,&parent,newEntities);
  else {
    for (size_t i = 0; i < newEntities.getSize(); i++) {
      int type = mesh->getType(newEntities[i]);
      if (type == 4) {
      	double val = apf::getScalar(field, parent, 0);
	apf::setScalar(field, newEntities[i], 0, val);
      }
    }
  }
}

void CavityTransfer::onCavity(
    EntityArray& oldElements,
    EntityArray& newEntities)
{
  apf::FieldShape *fs = apf::getShape(field);
  std::string name = fs->getName();
  bool fConstant = false;
  if (name == std::string("Constant_3"))
    fConstant = true;

  double sumOld = 0.;
  if (fConstant) {
    for (size_t i = 0; i < oldElements.getSize(); i++) {
      double val = apf::getScalar(field, oldElements[i], 0);
      sumOld += val*apf::measure(mesh, oldElements[i]);
      //sumOld += val;
    }
  }
  transfer(oldElements.getSize(),&(oldElements[0]),newEntities);

  double sumNew = 0.;
  double sumNewVal = 0.;
  if (fConstant) {
    for (size_t i = 0; i < newEntities.getSize(); i++) {
      int type = mesh->getType(newEntities[i]);
      if (type == 4) {
      	double val = apf::getScalar(field, newEntities[i], 0);
      	//double val = apf::measure(mesh, newEntities[i]);
      	sumNew += val*apf::measure(mesh, newEntities[i]);
      	//sumNew += val;
      }
    }
    for (size_t i = 0; i < newEntities.getSize(); i++) {
      int type = mesh->getType(newEntities[i]);
      if (type == 4) {
      	double val = apf::getScalar(field, newEntities[i], 0);
      	//double val = apf::measure(mesh, newEntities[i]);
	//apf::setScalar(field, newEntities[i], 0, val*sumOld/sumNew);
	if (sumNew == 0) continue;
	apf::setScalar(field, newEntities[i], 0, val*sumOld/sumNew);
      }
    }
    for (size_t i = 0; i < newEntities.getSize(); i++) {
      int type = mesh->getType(newEntities[i]);
      if (type == 4) {
      	double val = apf::getScalar(field, newEntities[i], 0);
      	//double val = apf::measure(mesh, newEntities[i]);
      	sumNewVal += val*apf::measure(mesh, newEntities[i]);
      }
    }
  }
}

/* hmm... could use multiple inheritance here, but that creates
   a "diamond problem":

        SolutionTransfer
          /           \
  LinearTranser    CavityTransfer
          \           /
        HighOrderTransfer

   and I don't want to reason about how
   stable it will be in this case.
   There are few transfer objects, so having duplicates is not too bad */
class HighOrderTransfer : public SolutionTransfer
{
  public:
    LinearTransfer verts;
    CavityTransfer others;
    HighOrderTransfer(apf::Field* f):
      verts(f),others(f)
    {
    }
    virtual bool hasNodesOn(int dimension)
    {
      return others.hasNodesOn(dimension);
    }
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert)
    {
      verts.onVertex(parent,xi,vert);
    }
    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities)
    {
      others.onRefine(parent,newEntities);
    }
    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities)
    {
      others.onCavity(oldElements,newEntities);
    }
};

SolutionTransfer* createFieldTransfer(apf::Field* f)
{
  apf::FieldShape* shape = apf::getShape(f);
  if (shape->hasNodesIn(0))
  {
    if (shape->getOrder() == 1)
      return new LinearTransfer(f);
    return new HighOrderTransfer(f);
  }
  /* turn this back on when IP fields are
     set up with shape functions */
  return new CavityTransfer(f);
}

SolutionTransfers::SolutionTransfers()
{
}

SolutionTransfers::~SolutionTransfers()
{
  for (size_t i = 0; i < transfers.size(); ++i)
    delete transfers[i];
}

void SolutionTransfers::add(SolutionTransfer* t)
{
  transfers.push_back(t);
}

bool SolutionTransfers::hasNodesOn(int dimension)
{
  for (size_t i = 0; i < transfers.size(); ++i)
    if (transfers[i]->hasNodesOn(dimension))
      return true;
  return false;
}

void SolutionTransfers::onVertex(
    apf::MeshElement* parent,
    Vector const& xi, 
    Entity* vert)
{
  for (size_t i = 0; i < transfers.size(); ++i)
    transfers[i]->onVertex(parent,xi,vert);
}

void SolutionTransfers::onRefine(
    Entity* parent,
    EntityArray& newEntities)
{
  for (size_t i = 0; i < transfers.size(); ++i)
    transfers[i]->onRefine(parent,newEntities);
}

void SolutionTransfers::onCavity(
    EntityArray& oldElements,
    EntityArray& newEntities)
{
  for (size_t i = 0; i < transfers.size(); ++i)
    transfers[i]->onCavity(oldElements,newEntities);
}

AutoSolutionTransfer::AutoSolutionTransfer(Mesh* m)
{
  for (int i = 0; i < m->countFields(); ++i)
  {
    apf::Field* f = m->getField(i);
    std::string name = f->getShape()->getName();
    if (name != std::string("Bezier"))
      this->add(createFieldTransfer(f));
  }
}

}
