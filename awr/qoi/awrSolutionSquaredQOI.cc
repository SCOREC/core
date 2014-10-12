/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrSolutionSquaredQOI.h"
#include <apfMesh.h>
#include <Teuchos_TestForException.hpp>

namespace awr {

/*****************************************************************************/
SolutionSquaredQOI::
SolutionSquaredQOI(apf::Mesh* m, const Teuchos::ParameterList& p) :
  QOI(m,p)
{
  validateParameters();
  /* assumes uniform mesh */
  init();
}

/*****************************************************************************/
void
SolutionSquaredQOI::
evaluateElementQOI(apf::MeshEntity* element,
                   apf::DynamicVector& f)
{
  f.setSize(num_nodes_);
}

/*****************************************************************************/
void
SolutionSquaredQOI::
init()
{
  apf::MeshIterator* elems = mesh_->begin(mesh_->getDimension());
  apf::MeshEntity* e = mesh_->iterate(elems);
  mesh_->end(elems);
  BasisUtils util(sol_,e,integration_order_);
  num_dims_ = util.getNumDims();
  num_nodes_ = util.getNumNodes();
  num_qp_ = util.getNumQP();
}

/*****************************************************************************/
void
SolutionSquaredQOI::
validateParameters()
{
  integration_order_ = params_.get("Integration Order",2); 
  std::string name = params_.get("Primal Solution Field Name","");
  sol_ = mesh_->findField(name.c_str());
  if (sol_ == NULL)
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: Poisson QOI: solution field with name \""
        << name << "\" does not exist\n");
}

/*****************************************************************************/
}
