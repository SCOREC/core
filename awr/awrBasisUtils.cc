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

BasisUtils::BasisUtils(apf::Mesh* m) :
  mesh_(m)
{
}

int BasisUtils::getNumDims()
{
  return mesh_->getDimension();
}

int BasisUtils::getNumQP(apf::MeshElement* elem, int order)
{
  return apf::countIntPoints(elem,order);
}

int BasisUtils::getNumNodes(apf::Field* f, apf::MeshEntity* elem)
{
  int type = mesh_->getType(elem);
  int num_nodes = f->getShape()->getEntityShape(type)->countNodes();
  return num_nodes;
}

void BasisUtils::getBF(apf::Element* elem, NodeQPScalar& bf)
{
}

void BasisUtils::getWBF(apf::Element* elem, NodeQPScalar& w_bf)
{
}

void BasisUtils::getGradBF(apf::Element* elem, NodeQPVector& grad_bf)
{
}

void BasisUtils::getWGradBF(apf::Element* elem, NodeQPVector& w_grad_bf)
{
}


void BasisUtils::
getFirstElem(apf::Field* f,
             apf::MeshEntity* elem,
             apf::MeshElement* me,
             apf::Element* fe)
{
  apf::MeshIterator* elems = mesh_->begin(mesh_->getDimension());
  elem = mesh_->iterate(elems);
  mesh_->end(elems);
  me = apf::createMeshElement(mesh_,elem);
  fe = apf::createElement(f,me);
}

}
