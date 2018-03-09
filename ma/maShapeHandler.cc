/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#include "maShapeHandler.h"
#include "maShape.h"
#include "maAdapt.h"
#include <apfShape.h>
#include <pcu_util.h>

namespace ma {

class LinearHandler : public ShapeHandler
{
  public:
    LinearHandler(Adapt* a)
    {
      mesh = a->mesh;
      sizeField = a->sizeField;
    }
    virtual double getQuality(Entity* e)
    {
      return measureElementQuality(mesh, sizeField, e);
    }
    virtual bool hasNodesOn(int dimension)
    {
      return dimension == 0;
    }
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert)
    {
      Vector point;
      apf::mapLocalToGlobal(parent,xi,point);
      mesh->setPoint(vert,0,point);
    }
  private:
    Mesh* mesh;
    SizeField* sizeField;
};

class QuadraticHandler : public ShapeHandler
{
  public:
    QuadraticHandler(Adapt* a)
    {
      mesh = a->mesh;
      st = createFieldTransfer(mesh->getCoordinateField());
    }
    ~QuadraticHandler()
    {
      delete st;
    }
    virtual double getQuality(Entity* e)
    {
      PCU_ALWAYS_ASSERT( mesh->getType(e) == apf::Mesh::TET );
      return measureQuadraticTetQuality(mesh,e);
    }
    virtual bool hasNodesOn(int dimension)
    {
      return st->hasNodesOn(dimension);
    }
    virtual void onVertex(
        apf::MeshElement* parent,
        Vector const& xi, 
        Entity* vert)
    {
      st->onVertex(parent,xi,vert);
    }
    virtual void onRefine(
        Entity* parent,
        EntityArray& newEntities)
    {
      st->onRefine(parent,newEntities);
    }
    virtual void onCavity(
        EntityArray& oldElements,
        EntityArray& newEntities)
    {
      st->onCavity(oldElements,newEntities);
    }
  private:
    Mesh* mesh;
    SolutionTransfer* st;
};

ShapeHandler* getShapeHandler(Adapt* a)
{
  apf::FieldShape* s = a->mesh->getShape();
  if (s->getOrder() == 1)
    return new LinearHandler(a);
  if (s->getOrder() == 2)
    return new QuadraticHandler(a);
  return 0;
}

}
