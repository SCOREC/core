/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <assert.h>
#include <PCU.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include "dwrUtils.h"
#include "dwrLinearSystem.h"
#include "dwrVectorL2QOI.h"
#include "dwrElasticityLHS.h"
#include "dwrElasticityProblem.h"

namespace dwr {

ElasticityProblem::ElasticityProblem()
{
  E = 0.0;
  nu = 0.0;
  quadratureDegree = 0;
  primal = NULL;
  dbc.setSize(0);
  owned_ = NULL;
  overlap_ = NULL;
}

ElasticityProblem::~ElasticityProblem()
{
  delete owned_;
  delete overlap_;
  delete ls_;
  apf::destroyGlobalNumbering(gn_);
}

void ElasticityProblem::validate()
{
  assert(E);
  assert(nu);
  assert(quadratureDegree);
  assert(primal);
  assert(dbc.getSize() > 0);
}

static Epetra_Map* createMap(int nc, apf::GlobalNumbering* n)
{
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(n,nodes);
  int nn = nodes.getSize();
  apf::DynamicArray<long long> dofIndices(nn*nc);
  for (int i=0; i < nn; ++i)
  {
    long global = apf::getNumber(n,nodes[i]);
    for (int j=0; j < nc; ++j)
      dofIndices[i*nc + j] = static_cast<long long>(global*nc + j);
  }
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  return new Epetra_Map(-1,dofIndices.getSize(),&dofIndices[0],0,comm);
}

static void number(
    apf::Field* f,
    long& nge,
    Epetra_Map** owned,
    Epetra_Map** overlap,
    apf::GlobalNumbering** gn)
{
  apf::Mesh* m = apf::getMesh(f);
  apf::FieldShape* fs = apf::getShape(f);
  apf::Numbering* n = apf::numberOwnedNodes(m,"dwr_node",fs);
  int nc = apf::countComponents(f);
  nge = static_cast<long>(apf::countNodes(n));
  nge *= static_cast<long>(nc);
  PCU_Add_Longs(&nge,1);
  *gn = apf::makeGlobal(n);
  *owned = createMap(nc,*gn);
  apf::synchronize(*gn);
  *overlap = createMap(nc,*gn);
}

void ElasticityProblem::setup()
{
  double t0 = PCU_Time();
  mesh_ = apf::getMesh(primal);
  number(primal,nge_,&owned_,&overlap_,&gn_);
  ls_ = new LinearSystem(nge_,owned_,overlap_);
  double t1 = PCU_Time();
  print("setup in %f seconds",t1-t0);
}

static void addKeToGlobalMatrix(
    LinearSystem* ls,
    apf::DynamicMatrix Ke,
    apf::MeshEntity* e,
    apf::Field* f,
    apf::GlobalNumbering* n)
{
  int nc = apf::countComponents(f);
  apf::NewArray<long> ngo;
  int nn = apf::getElementNumbers(n,e,ngo);
  for (int i=0; i < nn; ++i) {
  for (int c=0; c < nc; ++c) {
  for (int j=0; j < nn; ++j) {
  for (int d=0; d < nc; ++d) {
    int ii = i*nc + c;
    int jj = j*nc + d;
    long gi = ngo[i]*nc +c;
    long gj = ngo[j]*nc +d;
    ls->sumToMatrix(Ke(ii,jj),gj,gi);
  } } } }
}

static void addFeToGlobalVector(
    LinearSystem* ls,
    apf::DynamicVector Fe,
    apf::MeshEntity* e,
    apf::Field* f,
    apf::GlobalNumbering* n)
{
  int nc = apf::countComponents(f);
  apf::NewArray<long> ngo;
  int nn = apf::getElementNumbers(n,e,ngo);
  for (int i=0; i < nn; ++i) {
  for (int c=0; c < nc; ++c) {
    int ii = i*nc + c;
    long gi = ngo[i]*nc +c;
    ls->sumToVector(Fe(ii),gi);
  } }
}

static void applyBC(
    apf::Field* f,
    apf::DynamicArray<ElasticityDBC> const& dbc,
    apf::GlobalNumbering* gn,
    LinearSystem* ls)
{
  int nc = apf::countComponents(f);
  apf::Mesh* m = apf::getMesh(f);
  for (int i=0; i < dbc.getSize(); ++i)
  {
    int c = dbc[i].component;
    apf::ModelEntity* e = m->findModelEntity(dbc[i].dim,dbc[i].tag);
    apf::DynamicArray<apf::Node> nodes;
    apf::getNodesOnClosure(m,e,nodes);
    for (int n=0; n < nodes.getSize(); ++n)
    {
      long gid = apf::getNumber(gn,nodes[n])*nc + c;
      ls->diagonalizeMatrixRow(gid);
      ls->replaceToVector(0.0,gid);
    }
  }
}

void ElasticityProblem::assemble()
{
  double t0 = PCU_Time();
  VectorL2QOI qoi(quadratureDegree,primal);
  ElasticityLHS lhs(quadratureDegree,primal);
  lhs.setElasticModulus(E);
  lhs.setPoissonsRatio(nu);
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = mesh_->begin(mesh_->getDimension());
  while ((elem = mesh_->iterate(elems)))
  {
    apf::MeshElement* me = apf::createMeshElement(mesh_,elem);
    lhs.process(me);
    qoi.process(me);
    addKeToGlobalMatrix(ls_,lhs.Ke,elem,primal,gn_);
    addFeToGlobalVector(ls_,qoi.Fe,elem,primal,gn_);
    apf::destroyMeshElement(me);
  }
  mesh_->end(elems);
  applyBC(primal,dbc,gn_,ls_);
  ls_->completeMatrixFill();
  double t1 = PCU_Time();
  print("assembled in %f seconds",t1-t0);
}

static apf::Field* createDualField(apf::Field* p)
{
  apf::Mesh* m = apf::getMesh(p);
  const char* n = apf::getName(p);
  std::string name = std::string(n) + "_dual";
  int vt = apf::getValueType(p);
  apf::FieldShape* fs = apf::getShape(p);
  apf::Field* d = apf::createField(m,name.c_str(),vt,fs);
  return d;
}

static void attachSolution(
    apf::Mesh* m,
    apf::GlobalNumbering* gn,
    apf::Field* f,
    apf::DynamicVector const& sol)
{
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(gn,nodes);
  int nc = apf::countComponents(f);
  double v[nc];
  for (int i=0; i < nodes.getSize(); ++i)
  {
    for (int c=0; c < nc; ++c)
      v[c] = sol[i*nc + c];
    apf::setComponents(f,nodes[i].entity,nodes[i].node,v);
  }
  apf::synchronize(f);
}

void ElasticityProblem::solve()
{
  double t0 = PCU_Time();
  ls_->solve();
  apf::DynamicVector sol;
  ls_->getSolution(sol);
  dual = createDualField(primal);
  attachSolution(mesh_,gn_,dual,sol);
  double t1 = PCU_Time();
  print("solved in %f seconds",t1-t0);
}

apf::Field* ElasticityProblem::computeDual()
{
  print("solving linear elasticity dual problem");
  validate();
  setup();
  assemble();
  solve();
  return dual;
}

}
