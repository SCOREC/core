/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrBasisUtils.h"
#include "apfField.h"
#include "apfShape.h"

namespace awr {

BasisUtils::BasisUtils(apf::Mesh* m, apf::Field* f,
                       apf::MeshEntity* e, int o) :
  mesh_(m),
  sol_(f),
  elem_(e),
  order_(o)
{
  mesh_elem_ = apf::createMeshElement(mesh_,elem_);
  field_elem_ = apf::createElement(sol_,mesh_elem_);
}

BasisUtils::~BasisUtils()
{
  apf::destroyMeshElement(mesh_elem_);
  apf::destroyElement(field_elem_);
}

int BasisUtils::getNumDims()
{
  return mesh_->getDimension();
}

int BasisUtils::getNumQP()
{
  return apf::countIntPoints(mesh_elem_,order_);
}

int BasisUtils::getNumNodes()
{
  int type = mesh_->getType(elem_);
  int num_nodes = sol_->getShape()->getEntityShape(type)->countNodes();
  return num_nodes;
}

void BasisUtils::getBF(NodeQPScalar& bf)
{
}

void BasisUtils::getWBF(NodeQPScalar& w_bf)
{
}

void BasisUtils::getGradBF(NodeQPVector& grad_bf)
{
}

void BasisUtils::getWGradBF(NodeQPVector& w_grad_bf)
{
}

}
