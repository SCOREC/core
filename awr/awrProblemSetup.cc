/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrLinearSystem.h"
#include "awrQoI.h"
#include <PCU.h>
#include <apf.h>
#include <apfNumbering.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

namespace awr {

apf::Field* createAdjointField(apf::Mesh* m, apf::Field* p)
{
  const char* n = apf::getName(p);
  std::string name = std::string(n) + "_adj";
  int vt = apf::getValueType(p);
  apf::FieldShape* fs = apf::getShape(p);
  apf::Field* a = apf::createField(m,name.c_str(),vt,fs);
  return a;
}

apf::Numbering* createNumbering(apf::Mesh* m, apf::Field* a)
{
  apf::FieldShape* fs = apf::getShape(a);
  return apf::numberOwnedNodes(m,"node",fs);
}

long computeNumGlobalEqs(int nc, apf::Numbering* n)
{
  long nge = static_cast<long>(apf::countNodes(n));
  nge *= static_cast<long>(nc);
  PCU_Add_Longs(&nge,1);
  return nge;
}

apf::GlobalNumbering* globalizeNumbering(apf::Numbering* n)
{
  apf::GlobalNumbering* gn = apf::makeGlobal(n);
  return gn;
}

Epetra_Map* createMap(int nc, apf::GlobalNumbering* n)
{
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(n,nodes);
  int numNodes = nodes.getSize();
  apf::DynamicArray<long long> dofIndices(numNodes*nc);
  for (int i=0; i < numNodes; ++i)
  {
    long global = apf::getNumber(n,nodes[i]);
    std::cout << global << std::endl;
    for (int j=0; j < nc; ++j)
      dofIndices[i*nc + j] = global*nc + j;
  }
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  return new Epetra_Map(-1,dofIndices.getSize(),&dofIndices[0],0,comm);
}

void Problem::setup()
{
  double t0 = MPI_Wtime();
  validateProblemList(); /* pure virtual method */
  setPrimalField(); /* pure virtual method */
  adjoint_ = createAdjointField(mesh_,primal_);
  numComponents_ = apf::countComponents(adjoint_);
  numbering_ = createNumbering(mesh_,adjoint_);
  numGlobalEqs_ = computeNumGlobalEqs(numComponents_,numbering_);
  globalNumbering_ = globalizeNumbering(numbering_);
  Epetra_Map* owned = createMap(numComponents_,globalNumbering_);
  apf::synchronize(globalNumbering_);
  Epetra_Map* overlap = createMap(numComponents_,globalNumbering_);
  qoi_ = createQoI(qoiList_,mesh_,primal_);
  ls_ = new LinearSystem(numGlobalEqs_,owned,overlap);
  double t1 = MPI_Wtime();
  print("set up in %f seconds",t1-t0);
}

}
