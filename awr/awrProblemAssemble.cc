/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrQoI.h"
#include "awrLinearSystem.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfDynamicMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

void addKeToGlobalMatrix(
    LinearSystem* ls,
    apf::DynamicMatrix Ke,
    apf::MeshEntity* e,
    apf::Field* f,
    apf::GlobalNumbering* n)
{
  int nc = apf::countComponents(f);
  apf::NewArray<long> ngo;
  int nn = apf::getElementNumbers(n,e,ngo);
  int nd = nn*nc;
  for (int i=0; i < nd; ++i)
  {
    for (int j=0; j < nd; ++j)
    {
      long gi = ngo[i];
      long gj = ngo[j];
      ls->sumToMatrix(Ke(i,j),gi,gj);
    }
  }
}

void addFeToGlobalVector(
    LinearSystem* ls,
    apf::DynamicVector Fe,
    apf::MeshEntity* e,
    apf::Field* f,
    apf::GlobalNumbering* n)
{
  int nc = apf::countComponents(f);
  apf::NewArray<long> ngo;
  int nn = apf::getElementNumbers(n,e,ngo);
  int nd = nn*nc;
  for (int i=0; i < nd; ++i)
  {
    long gi = ngo[i];
    ls->sumToVector(Fe(i),gi);
  }
}

void Problem::assemble()
{
  createIntegrator();
  qoi_->createIntegrator();
  apf::MeshEntity* elem;
  apf::DynamicMatrix Ke;
  apf::DynamicVector Fe;
  apf::MeshIterator* elements =
    mesh_->begin(mesh_->getDimension());
  while ((elem = mesh_->iterate(elements)))
  {
    processKe(elem,Ke);
    qoi_->processFe(elem,Fe);
    addKeToGlobalMatrix(ls_,Ke,elem,adjoint_,globalNumbering_);
    addFeToGlobalVector(ls_,Fe,elem,adjoint_,globalNumbering_);
  }
  mesh_->end(elements);
  ls_->completeMatrixFill();
}

}
