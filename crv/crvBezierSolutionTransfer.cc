/*
1;3409;0c * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apfElement.h>

#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"
#include "crvQuality.h"
#include "crvShape.h"
#include "crvTables.h"
#include <maShapeHandler.h>
#include <maSolutionTransfer.h>
#include <maShape.h>
#include <mth_def.h>
#include <math.h>
#include <pcu_util.h>

namespace crv {
class LinearTransfer;
class CavityTransfer;

CrvBezierFieldTransfer::~CrvBezierFieldTransfer()
{
}

void CrvBezierFieldTransfer::onVertex(
        apf::MeshElement*,
        Vector const&, 
        Entity*)
{
}

void CrvBezierFieldTransfer::onRefine(
    Entity*,
    EntityArray&)
{
}

void CrvBezierFieldTransfer::onCavity(
    EntityArray&,
    EntityArray&)
{
}

class HigherOrderTransfer : public CrvBezierField
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
    virtual void convertToBezierFields(
        int dimension,
        EntityArray& newEntities)
    {
      //apf::FieldShape *fs = others.shape;
      //std::string name = fs->getName();
      //if (name == std::string("Bezier")) {
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
};

CrvBezierFieldTransfer::CrvBezierFieldTransfer(Mesh* m)
{
  for (int i = 0; i < m->countFields(); ++i)
  {
    apf::Field* f = m->getField(i);
    printf("Fields considered, %s\n", apf::getName(f));
    this->add(createFieldTransfer(f));
  }
}


}

