/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrLinearSystem.h"
#include <PCU.h>
#include <apf.h>
#include <apfNumbering.h>

namespace awr {

void Problem::createAdjointField()
{
  const char* n = apf::getName(primal_);
  std::string name = std::string(n) + "_adj";
  int vt = apf::getValueType(primal_);
  apf::FieldShape* fs = apf::getShape(primal_);
  adjoint_ = apf::createField(mesh_,name.c_str(),vt,fs);
  numComponents_ = apf::countComponents(adjoint_);
  /*********************************************/
  apf::MeshEntity* v;
  apf::MeshIterator* vertices = mesh_->begin(0);
  while ((v = mesh_->iterate(vertices)))
  {
    apf::setScalar(adjoint_,v,0,1.0);
  }
  mesh_->end(vertices);
  /*********************************************/
}

void Problem::createNumbering()
{
  apf::FieldShape* fs = apf::getShape(primal_);
  numbering_ = apf::numberOwnedNodes(mesh_,"dof",fs);
}

void Problem::computeNumGlobalEqs()
{
  numGlobalEqs_ = static_cast<long>(apf::countNodes(numbering_));
  numGlobalEqs_ *= static_cast<long>(numComponents_);
  PCU_Add_Longs(&numGlobalEqs_,1);
}

void Problem::globalizeNumbering()
{
  globalNumbering_ = apf::makeGlobal(numbering_);
  apf::synchronize(globalNumbering_);
}

void Problem::setup()
{
  /* pure virtual method */
  validateProblemList();
  /* pure virtual method */
  setPrimalField();
  createAdjointField();
  createNumbering();
  computeNumGlobalEqs();
  globalizeNumbering();
  ls_ = new LinearSystem(numGlobalEqs_);
}

}
