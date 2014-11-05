/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRLINEARSYSTEM_H
#define AWRLINEARSYSTEM_H

#include "awrTrilinosTypes.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace awr {

class LinearSystem
{
  public:
    LinearSystem(GO n,const Teuchos::RCP<Teuchos::ParameterList> p);
    virtual ~LinearSystem() {};
    GO getNumGlobalUnknowns() { return numGlobalUnknowns_; }
    LO getNumLocalUnknowns() { return numLocalUnknowns_; }
    GO mapLocalToGlobal(LO i);
    void insertToVector(double v, GO i);
    void insertToMatrix(double v, GO i, GO j);
    void sumToVector(double v, GO i);
    void sumToMatrix(double v, GO i, GO j);
    void completeMatrixFill();
    void describeMatrix();
    void describeVector();
    void solve();
  private:
    GO numGlobalUnknowns_;
    LO numLocalUnknowns_;
    Teuchos::RCP<MatrixT> A_;
    Teuchos::RCP<VectorT> b_;
    Teuchos::RCP<const MapT> map_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    LinearSystem(const LinearSystem&);
    LinearSystem& operator=(const LinearSystem&);
};

}

#endif
