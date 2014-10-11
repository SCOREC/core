/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrPoissonLHS.h"
#include <apfMesh.h>
#include <Teuchos_TestForException.hpp>

namespace awr {

/*****************************************************************************/
PoissonLHS::
PoissonLHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  LHS(m,p)
{
  validateParameters();
  /* assumes uniform mesh */
  init();
}

/*****************************************************************************/
void
PoissonLHS::
evaluateElementLHS(apf::MeshEntity* element,
                   apf::DynamicMatrix& k)
{
  k.setSize(num_nodes_,num_nodes_);
  for (int a=0; a < k.getRows(); ++a)
  {
    for (int b=0; b < k.getColumns(); ++b)
    {
      k(a,b) = 0.0;
      for (int qp=0; qp < num_qp_; ++qp)
      for (int i=0; i < num_dims_; ++i)
        k(a,b) += grad_bf_[a][qp][i] * w_grad_bf_[b][qp][i];
    }
  }
}

/*****************************************************************************/
void
PoissonLHS::
init()
{
  apf::MeshIterator* elems = mesh_->begin(mesh_->getDimension());
  apf::MeshEntity* e = mesh_->iterate(elems);
  mesh_->end(elems);
  BasisUtils util(sol_,e,integration_order_);
  num_dims_ = util.getNumDims();
  num_nodes_ = util.getNumNodes();
  num_qp_ = util.getNumQP();
  util.getGradBF(grad_bf_);
  util.getWGradBF(w_grad_bf_);
}

/*****************************************************************************/
void
PoissonLHS::
validateParameters()
{
  integration_order_ = params_.get("Integration Order",2); 
  std::string name = params_.get("Primal Solution Field Name","");
  sol_ = mesh_->findField(name.c_str());
  if (sol_ == NULL)
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: Poisson LHS: solution field with name \""
        << name << "\" does not exist\n");
}

/*****************************************************************************/
}
