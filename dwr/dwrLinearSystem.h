/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_LINEARSYSTEM_H
#define DWR_LINEARSYSTEM_H

class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_MultiVector;

namespace apf {
class DynamicVector;
}

namespace dwr {

class LinearSystem
{
  private:
    typedef int LO;
    typedef long long GO;
  public:
    LinearSystem(GO n, Epetra_Map* owned, Epetra_Map* overlap);
    ~LinearSystem();
    void sumToVector(double v, GO i);
    void replaceToVector(double v, GO i);
    void sumToMatrix(double v, GO i, GO j);
    void diagonalizeMatrixRow(GO i);
    void completeMatrixFill();
    void solve();
    void getSolution(apf::DynamicVector& sol);
  private:
    GO numGlobalEqs_;
    Epetra_Map* ownedMap_;
    Epetra_Map* overlapMap_;
    Epetra_CrsMatrix* A_;
    Epetra_MultiVector* x_;
    Epetra_MultiVector* b_;
};

}

#endif
