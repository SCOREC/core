/*
1;3409;0c * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apfElement.h>

#include "crv.h"
#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"
#include "crvQuality.h"
#include "crvShape.h"
#include "crvTables.h"
#include <maMesh.h>
#include <maShapeHandler.h>
#include <maSolutionTransfer.h>
#include <maShape.h>
#include <mth_def.h>
#include <math.h>
#include <pcu_util.h>

namespace crv {

class CrvBezierSolutionTransfer : public ma::SolutionTransfer
{
  public:
    ma::LinearTransfer verts;
    ma::CavityTransfer others;
    CrvBezierSolutionTransfer(apf::Field* f):
      verts(f),others(f)
    {
    }
    virtual bool hasNodesOn(int dimension)
    {
      return others.hasNodesOn(dimension);
    }
    virtual void onVertex(
        apf::MeshElement* parent,
        ma::Vector const& xi,
        ma::Entity* vert)
    {
      verts.onVertex(parent,xi,vert);
    }
    virtual void onRefine(
        ma::Entity* parent,
        ma::EntityArray& newEntities)
    {
      others.onRefine(parent,newEntities);
      apf::FieldShape *shape = others.shape;
      apf::Mesh *m = others.mesh;
      int type = m->getType(newEntities[0]);
      int td = apf::Mesh::typeDimension[type];
      int order = shape->getOrder();
      int n = shape->getEntityShape(
          apf::Mesh::simplexTypes[td])->countNodes();
      int ne = shape->countNodesOn(
          apf::Mesh::simplexTypes[td]);
      apf::NewArray<double> c;
      crv::getBezierTransformationCoefficients(order,
          apf::Mesh::simplexTypes[td], c);

      for( size_t i = 0; i < newEntities.getSize(); i++) {
        crv::convertInterpolationFieldPoints(newEntities[i],
            others.field, n, ne, c);
      }
    }
    virtual void onCavity(
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {
      others.onCavity(oldElements,newEntities);
    }
    /*
  protected:
    void convertToBezierFields(
        int dimension,
        EntityArray& newEntities)
    {
      //apf::FieldShape *fs = others.shape;
      //std::string name = fs->getName();
      //if (name == std::string("Bezier")) {
      int td = apf::Mesh::typeDimension[]
      int order = shape->getOrder();
      int n = shape->getEntityShape(
          apf::Mesh::simplexTypes[dimension])->countNodes();
      int ne = shape->countNodesOn(
          apf::Mesh::simplexTypes[dimension]);
      apf::NewArray<double> c;
      crv::getBezierTransformationCoefficients(order,
          apf::Mesh::simplexTypes[dimension], c);

      for( size_t i = 0; i < newEntities.getSize(); i++) {
        crv::convertInterpolationFieldPoints(newEntities[i],
            field, n, ne, c);
      }
    }
    */
};

static ma::SolutionTransfer* createBezierSolutionTransfer(apf::Field* f)
{
  return new CrvBezierSolutionTransfer(f);
}

ma::SolutionTransfer* setBezierSolutionTransfers(
    const std::vector<apf::Field*>& fields)
{
  ma::SolutionTransfers* st = new ma::SolutionTransfers();
  for (std::size_t i = 0; i < fields.size(); i++) {
    st->add(createBezierSolutionTransfer(fields[i]));
  }
  return st;
}

}
