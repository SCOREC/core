/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrPoissonRHS.h"
#include "Teuchos_TestForException.hpp"

namespace awr {

/*****************************************************************************/
PoissonRHS::
PoissonRHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  RHS(m,p)
{
  validateParameters();
  /* assumes uniform mesh */
  init();
}

/*****************************************************************************/
void
PoissonRHS::
evaluateElementRHS(apf::MeshEntity* element,
                   apf::DynamicMatrix& k)
{
}

/*****************************************************************************/
void
PoissonRHS::
init()
{
  BasisUtils util(mesh_);
  num_dims_ = util.getNumDims();
  apf::MeshIterator* elems = mesh_->begin(num_dims_);
  apf::MeshEntity* e = mesh_->iterate(elems);
  mesh_->end(elems);
  num_nodes_ = util.getNumNodes(sol_,e);
  apf::MeshElement* me = apf::createMeshElement(mesh_,e);
  num_qp_ = util.getNumQP(me,integration_order_);
  apf::Element* fe = createElement(sol_,me);
  util.getGradBF(fe,grad_bf_);
  util.getWGradBF(fe,w_grad_bf_);
  apf::destroyMeshElement(me);
  apf::destroyElement(fe);
}

/*****************************************************************************/
void
PoissonRHS::
validateParameters()
{
  integration_order_ = params_.get("Integration Order",2); 
  std::string name = params_.get("Primal Solution Field Name","");
  sol_ = mesh_->findField(name.c_str());
  if (sol_ == NULL)
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: Poisson RHS: solution field with name \""
        << name << "\" does not exist\n");
                             
}

/*****************************************************************************/
}
