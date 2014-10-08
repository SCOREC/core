/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrBasisUtils.h"
#include "apfField.h"
#include "apfShape.h"
#include "apfVectorElement.h"

namespace awr {

BasisUtils::BasisUtils(apf::Field* f, apf::MeshEntity* e, int o) :
  sol_(f),
  elem_(e),
  order_(o)
{
  mesh_ = getMesh(sol_);
  mesh_elem_ = apf::createMeshElement(mesh_,elem_);
  field_elem_ = apf::createElement(sol_,mesh_elem_);
  int type = mesh_->getType(elem_);
  num_nodes_ = sol_->getShape()->getEntityShape(type)->countNodes();
  num_dims_ = mesh_->getDimension();
  num_qp_ = apf::countIntPoints(mesh_elem_,order_);
}

BasisUtils::~BasisUtils()
{
  apf::destroyMeshElement(mesh_elem_);
  apf::destroyElement(field_elem_);
}

void BasisUtils::getBF(NodeQPScalar& bf)
{
  bf.setSize(num_nodes_);
  for (int n=0; n < num_nodes_; ++n)
    bf[n].setSize(num_qp_);
  for (int qp=0; qp < num_qp_; ++qp)
  {
    apf::Vector3 param;
    apf::getIntPoint(mesh_elem_,order_,qp,param);
    int type = mesh_->getType(elem_);
    apf::NewArray<double> values;
    sol_->getShape()->getEntityShape(type)->getValues(param,values);
    for (int n=0; n < num_nodes_; ++n)
      bf[n][qp] = values[n];
  }
}

void BasisUtils::getWBF(NodeQPScalar& w_bf)
{
  this->getBF(w_bf);
  for (int qp=0; qp < num_qp_; ++qp)
  {
    apf::Vector3 param;
    apf::getIntPoint(mesh_elem_,order_,qp,param);
    double w = apf::getIntWeight(mesh_elem_,order_,qp);
    double j = apf::getDV(mesh_elem_,param);
    for (int n=0; n < num_qp_; ++n)
      w_bf[n][qp] *= w*j;
  }
}

void BasisUtils::getGradBF(NodeQPVector& grad_bf)
{
  grad_bf.setSize(num_nodes_);
  for (int n=0; n < num_nodes_; ++n)
    grad_bf[n].setSize(num_nodes_);
  for (int n=0; n < num_nodes_; ++n)
  for (int qp=0; qp < num_qp_; ++qp)
    grad_bf[n][qp].setSize(num_dims_);
  for (int qp=0; qp < num_qp_; ++qp)
  {
    apf::Vector3 param;
    apf::getIntPoint(mesh_elem_,order_,qp,param);
    apf::Matrix3x3 j;
    mesh_elem_->getJacobian(param,j);
    apf::Matrix3x3 j_inv = invert(j);
    int type = mesh_->getType(elem_);
    apf::NewArray<apf::Vector3> grads;
    sol_->getShape()->getEntityShape(type)->
      getLocalGradients(param,grads);
    for (int n=0; n < num_nodes_; ++n)
    {
      apf::Vector3 global_grad;
      global_grad = j_inv*grads[n];
      for (int i=0; i < num_dims_; ++i)
        grad_bf[n][qp][i] = global_grad[i];
    }
  }
}

void BasisUtils::getWGradBF(NodeQPVector& w_grad_bf)
{
  this->getGradBF(w_grad_bf);
  for (int qp=0; qp < num_qp_; ++qp)
  {
    apf::Vector3 param;
    apf::getIntPoint(mesh_elem_,order_,qp,param);
    double w = apf::getIntWeight(mesh_elem_,order_,qp);
    double j = apf::getDV(mesh_elem_,param);
    for (int n=0; n < num_qp_; ++n)
    for (int i=0; i < num_dims_; ++i)
      w_grad_bf[n][qp][i] *= w*j;
  }
}

}
