/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apfElement.h>

#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvQuality.h"
#include <apfShape.h>
#include <apfNumbering.h>
#include <maAffine.h>
#include <maMap.h>
#include <maShapeHandler.h>
#include <maSolutionTransfer.h>
#include <cassert>
#include <float.h>

namespace crv {

class FieldTransfer : public ma::SolutionTransfer
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

class BezierTransfer : public FieldTransfer
{
  public:
    BezierTransfer(apf::Field* f):
      FieldTransfer(f), minDim(0)
    {
    }
    void transferToNodeIn(
        apf::Element* elem,
        apf::Node const& node,
        ma::Vector const& elemXi)
    {
      apf::getComponents(elem,elemXi,&(value[0]));
      apf::setComponents(field,node.entity,node.node,&(value[0]));
    }
    int getBestElement(
        int n,
        apf::Element** elems,
        ma::Affine* elemInvMaps,
        ma::Vector const& point,
        ma::Vector& bestXi)
    {
      double bestValue = -DBL_MAX;
      int bestI = 0;
      for (int i = 0; i < n; ++i)
      {
        ma::Vector xi = elemInvMaps[i] * point;
        double value = ma::getInsideness(mesh,apf::getMeshEntity(elems[i]),xi);
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
        int n,
        apf::Element** elems,
        ma::Affine* elemInvMaps,
        apf::Node const& node)
    {
      ma::Vector xi;
      shape->getNodeXi(mesh->getType(node.entity),node.node,xi);
      ma::Affine childMap = ma::getMap(mesh,node.entity);
      ma::Vector point = childMap * xi;
      ma::Vector elemXi;
      int i = getBestElement(n,elems,elemInvMaps,point,elemXi);
      transferToNodeIn(elems[i],node,elemXi);
    }
    void transfer(
        int n,
        ma::Entity** cavity,
        ma::EntityArray& newEntities)
    {
      if (getDimension(mesh, cavity[0]) < minDim)
        return;
      apf::NewArray<apf::Element*> elems(n);
      for (int i = 0; i < n; ++i)
        elems[i] = apf::createElement(field,cavity[i]);
      apf::NewArray<ma::Affine> elemInvMaps(n);
      for (int i = 0; i < n; ++i)
        elemInvMaps[i] = invert(ma::getMap(mesh,cavity[i]));
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
    virtual void onRefine(
        ma::Entity* parent,
        ma::EntityArray& newEntities)
    {
      transfer(1,&parent,newEntities);
    }
    virtual void onCavity(
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {
      transfer(oldElements.getSize(),&(oldElements[0]),newEntities);
    }
  private:
    int minDim;
};

class BezierHandler : public ma::ShapeHandler
{
  public:
    BezierHandler(ma::Mesh* m)
    {
      mesh = m;
      bt = new BezierTransfer(mesh->getCoordinateField());
    }
    ~BezierHandler()
    {
      delete bt;
    }
    virtual double getQuality(apf::MeshEntity* e)
    {
      assert( mesh->getType(e) == apf::Mesh::TRIANGLE ||
          mesh->getType(e) == apf::Mesh::TET);
      return crv::getQuality(mesh,e);
    }
    virtual bool hasNodesOn(int dimension)
    {
      return bt->hasNodesOn(dimension);
    }
    virtual void onVertex(
        apf::MeshElement* parent,
        apf::Vector3 const& xi,
        apf::MeshEntity* vert)
    {
      bt->onVertex(parent,xi,vert);
    }
    virtual void onRefine(
        apf::MeshEntity* parent,
        ma::EntityArray& newEntities)
    {
      bt->onRefine(parent,newEntities);
    }
    virtual void onCavity(
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {
      bt->onCavity(oldElements,newEntities);
    }
  private:
    ma::Mesh* mesh;
    BezierTransfer* bt;
};

ma::ShapeHandler* getShapeHandler(ma::Mesh* m)
{
  return new BezierHandler(m);
}

}
