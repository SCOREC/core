/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "dwrLinearSystem.h"
#include <AztecOO.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <apfDynamicVector.h>

namespace dwr {

LinearSystem::LinearSystem(GO n, Epetra_Map* owned, Epetra_Map* overlap) :
  numGlobalEqs_(n),
  ownedMap_(owned),
  overlapMap_(overlap)
{
  A_ = new Epetra_CrsMatrix(Copy,*overlapMap_,n);
  x_ = new Epetra_MultiVector(*ownedMap_,/*num vectors=*/1);
  b_ = new Epetra_MultiVector(*overlapMap_,/*num vectors=*/1);
}

LinearSystem::~LinearSystem()
{
  delete A_;
  delete x_;
  delete b_;
}

void LinearSystem::sumToVector(double v, GO i)
{
  int err = b_->SumIntoGlobalValue(i,/*vec idx=*/0,v);
  assert(err == 0);
}

void LinearSystem::replaceToVector(double v, GO i)
{
  b_->ReplaceGlobalValue(i,/*vec idx=*/0,v);
}

void LinearSystem::sumToMatrix(double v, GO i, GO j)
{
  double val[1]; val[0] = v;
  GO col[1]; col[0] = j;
  int err = A_->SumIntoGlobalValues(i,1,val,col);
  if (err != 0) {
    err = A_->InsertGlobalValues(i,1,val,col);
    assert(err == 0);
  }
}

void LinearSystem::diagonalizeMatrixRow(GO i)
{
  int n;
  double* v;
  GO* j;
  A_->ExtractGlobalRowView(i,n,v,j);
  for (int k=0; k < n; ++k)
    v[k] = 0.0;
  this->sumToMatrix(1.0,i,i);
}

void LinearSystem::completeMatrixFill()
{
  A_->FillComplete();
}

void LinearSystem::getSolution(apf::DynamicVector& sol)
{
  Epetra_Import importer(*overlapMap_,*ownedMap_);
  Epetra_MultiVector x(*overlapMap_,/*num vectors=*/1);
  assert(x.Import(*x_,importer,Add) == 0);
  double** s;
  x.ExtractView(&s);
  sol.setSize(x.MyLength());
  for (int i=0; i < x.MyLength(); ++i)
    sol[i] = s[0][i];
}

void LinearSystem::solve()
{
  Epetra_Export exporter(*overlapMap_,*ownedMap_);
  Epetra_CrsMatrix A(Copy,*ownedMap_,numGlobalEqs_);
  Epetra_MultiVector b(*ownedMap_,/*num vectors=*/1);
  assert(A.Export(*A_,exporter,Add) == 0);
  assert(b.Export(*b_,exporter,Add) == 0);
  A.FillComplete();
  Epetra_LinearProblem problem(&A,x_,&b);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond,AZ_Jacobi);
  solver.SetAztecOption(AZ_output,AZ_none);
  solver.Iterate(1000,1.0e-8);
}

}
