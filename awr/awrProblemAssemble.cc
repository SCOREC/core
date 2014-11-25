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

void rejectBCList(const char* msg)
{
  print("bc parameter %s\n",msg);
  abort();
}

void parseBCList(
    ParameterList& p,
    Teuchos::Array<int>& tags,
    Teuchos::Array<int>& dims)
{
  if (! p.isParameter("Geometric Tags"))
    rejectBCList("Geometric Tags not set");
  if (! p.isParameter("Geometric Dims"))
    rejectBCList("Geometric Dims not set");
  tags = p.get<Teuchos::Array<int> >("Geometric Tags");
  dims = p.get<Teuchos::Array<int> >("Geometric Dims");
  assert(tags.size() == dims.size());
}

void addBCToLinearSystem(
    apf::Mesh* m,
    LinearSystem* ls,
    apf::GlobalNumbering* gn,
    int dim,
    int tag)
{
  apf::ModelEntity* e = m->findModelEntity(dim,tag);
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodesOnClosure(m,e,nodes);
  for (int i=0; i < nodes.getSize(); ++i)
  {
    long gid = apf::getNumber(gn,nodes[i]);
    ls->diagonalizeMatrixRow(gid);
    ls->replaceToVector(0.0,gid);
  }
}

void applyBC(
    apf::Mesh* m,
    ParameterList& p,
    apf::GlobalNumbering* gn,
    LinearSystem* ls)
{
  Teuchos::Array<int> dims;
  Teuchos::Array<int> tags;
  parseBCList(p,tags,dims);
  for (int i=0; i < tags.size(); ++i)
    addBCToLinearSystem(m,ls,gn,dims[i],tags[i]);
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
  applyBC(mesh_,bcList_,globalNumbering_,ls_);
  ls_->completeMatrixFill();
}

}
