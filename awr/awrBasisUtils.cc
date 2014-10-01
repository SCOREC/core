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
  int num_nodes = this->getNumNodes();
  int num_qp = this->getNumQP();
  bf.setSize(num_nodes);
  for (int n=0; n < num_nodes; ++n)
    bf[n].setSize(num_qp);
  apf::NewArray<double> values;
  apf::Vector3 param;
  for (int qp=0; qp < num_qp; ++qp)
  {
    apf::getIntPoint(mesh_elem_,order_,qp,param);
    int type = mesh_->getType(elem_);
    sol_->getShape()->getEntityShape(type)->getValues(param,values);
    for (int n=0; n < num_nodes; ++n)
      bf[n][qp] = values[n];
  }
}

void BasisUtils::getWBF(NodeQPScalar& w_bf)
{
}

void BasisUtils::getGradBF(NodeQPVector& grad_bf)
{
  int num_nodes = this->getNumNodes();
  int num_qp = this->getNumQP();
  int num_dims = this->getNumDims();
  grad_bf.setSize(num_nodes);
  for (int n=0; n < num_nodes; ++n)
    grad_bf[n].setSize(num_nodes);
  for (int n=0; n < num_nodes; ++n)
  for (int qp=0; qp < num_qp; ++qp)
    grad_bf[n][qp].setSize(num_dims);
  apf::NewArray<apf::Vector3> grads;
  apf::Vector3 param;
  for (int qp=0; qp < num_qp; ++qp)
  {
    apf::getIntPoint(mesh_elem_,order_,qp,param);
    int type = mesh_->getType(elem_);
    sol_->getShape()->getEntityShape(type)->
      getLocalGradients(param,grads);
    for (int n=0; n < num_nodes; ++n)
    for (int i=0; i < num_dims; ++i)
      grad_bf[n][qp][i] = grads[n][i];
  }
}

void BasisUtils::getWGradBF(NodeQPVector& w_grad_bf)
{
}

}
