/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSolutionTransfer.h"
#include "maSolutionTransferHelper.h"
#include <apfShape.h>
#include <apfNumbering.h>
#include <float.h>

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

class FieldTransfer : public SolutionTransfer
{
  public:
    FieldTransfer(apf::Field* f)
    {
      field = f;
      mesh = apf::getMesh(f);
      shape = apf::getShape(f);
      value.allocate(apf::countComponents(f));
    }
    /* hmm... in vs. on ... probably the ma:: signature
       should change, it has the least users */
    virtual bool hasNodesOn(int dimension)
    {
      return shape->hasNodesIn(dimension);
    }
    apf::Field* field;
    apf::Mesh* mesh;
    apf::FieldShape* shape;
    apf::NewArray<double> value;
};

class LinearTransfer : public FieldTransfer
{
  public:
    LinearTransfer(apf::Field* f):
      FieldTransfer(f)
    {
    }
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert)
    {
      apf::Element* e = apf::createElement(field,parent);
      apf::getComponents(e,xi,&(value[0]));
      apf::setComponents(field,vert,0,&(value[0]));
      apf::destroyElement(e);
    }
};

class CavityTransfer : public FieldTransfer
{
  public:
    CavityTransfer(apf::Field* f):
      FieldTransfer(f)
    {
      minDim = getMinimumDimension(getShape(f));
    }
    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities)
    {
      transfer(field, &(value[0]), 1,&parent,newEntities);
    }
    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities)
    {
      transfer(field, &(value[0]),
           oldElements.getSize(),&(oldElements[0]),newEntities);
    }
  private:
    int minDim;
};

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
    this->add(createFieldTransfer(f));
  }
}

}
