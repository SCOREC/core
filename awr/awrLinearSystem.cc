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
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>

namespace awr
{

LinearSystem::LinearSystem(GO n, Epetra_Map* map) :
  numGlobalEqs_(n),
  map_(map)
{
  numLocalEqs_ = map_->NumMyElements();
  A_ = new Epetra_CrsMatrix(Copy,*map_,n);
  x_ = new Epetra_MultiVector(*map_,/*num vectors=*/1);
  b_ = new Epetra_MultiVector(*map_,/*num vectors=*/1);
}

LinearSystem::~LinearSystem()
{
  delete map_;
  delete A_;
  delete x_;
  delete b_;
}

GO LinearSystem::mapLIDtoGID(LO lid)
{
  GO* myGlobalElements = NULL;
  myGlobalElements = map_->MyGlobalElements64();
  return myGlobalElements[lid];
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
  A_->InsertGlobalValues(i,1,val,col);
}

void LinearSystem::diagonalizeMatrixRow(GO i)
{
  int n;
  double* v;
  GO* j;
  A_->ExtractGlobalRowView(i,n,v,j);
  for (int i=0; i < n; ++i)
    v[i] = 0.0;
  this->sumToMatrix(1.0,i,i);
}

void LinearSystem::completeMatrixFill()
{
  A_->FillComplete();
}

void LinearSystem::solve()
{
  Epetra_LinearProblem problem(A_,x_,b_);
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond,AZ_Jacobi);
  solver.Iterate(1000,1.0e-8);
}

}
