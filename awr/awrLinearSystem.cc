/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrLinearSystem.h"
#include <AztecOO.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>

namespace awr
{

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
  delete ownedMap_;
  delete overlapMap_;
}

void LinearSystem::sumToVector(double v, GO i)
{
  b_->SumIntoGlobalValue(i,/*vec idx=*/0,v);
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
  if (err != 0)
    A_->InsertGlobalValues(i,1,val,col);
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

double* LinearSystem::getSolution()
{
  double** sol;
  x_->ExtractView(&sol);
  return sol[0];
}

void LinearSystem::solve()
{
  Epetra_Import importer(*ownedMap_,*overlapMap_);
  Epetra_CrsMatrix A(Copy,*ownedMap_,numGlobalEqs_);
  Epetra_MultiVector b(*ownedMap_,/*num vectors=*/1);
  assert(A.Import(*A_,importer,Insert) == 0);
  assert(b.Import(*b_,importer,Insert) == 0);

  /* why isn't import add working?
     incorrect owned map? */
  b_->Print(std::cout);
  b.Print(std::cout);

  Epetra_LinearProblem problem(&A,x_,&b);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond,AZ_Jacobi);
  solver.SetAztecOption(AZ_output,AZ_none);
  solver.Iterate(1000,1.0e-8);
}

}
