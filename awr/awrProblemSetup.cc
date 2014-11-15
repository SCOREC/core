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
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

namespace awr {

void Problem::createAdjointField()
{
  const char* n = apf::getName(primal_);
  std::string name = std::string(n) + "_adj";
  int vt = apf::getValueType(primal_);
  apf::FieldShape* fs = apf::getShape(primal_);
  adjoint_ = apf::createField(mesh_,name.c_str(),vt,fs);
  numComponents_ = apf::countComponents(adjoint_);
  /* this is just so apf::writeVTKFiles doesn't blow
   * up for now */
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
  numbering_ = apf::numberOwnedNodes(mesh_,"node",fs);
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

Epetra_Map* createMap(int nc, apf::GlobalNumbering* numbering)
{
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(numbering,nodes);
  int numOverlapNodes = nodes.getSize();
  apf::DynamicArray<long long> dofIndices(numOverlapNodes*nc);
  for (int i=0; i < numOverlapNodes; ++i)
  {
    long global = apf::getNumber(numbering,nodes[i]);
    for (int j=0; j < nc; ++j)
      dofIndices[i*nc + j] = global*nc + j;
  }
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  return new Epetra_Map(-1,dofIndices.getSize(),&dofIndices[0],0,comm);
}

void Problem::setup()
{
  validateProblemList(); /* pure virtual method */
  setPrimalField(); /* pure virtual method */
  createAdjointField();
  createNumbering();
  computeNumGlobalEqs();
  globalizeNumbering();
  ls_ = new LinearSystem(numGlobalEqs_,
      createMap(numComponents_,globalNumbering_));
}

}
