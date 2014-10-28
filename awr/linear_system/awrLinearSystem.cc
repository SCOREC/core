/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrLinearSystem.h"
#include <Tpetra_DefaultPlatform.hpp>

namespace awr {

LinearSystem::
LinearSystem(GO n, const Teuchos::RCP<Teuchos::ParameterList> p) :
  numGlobalUnknowns_(n),
  params_(p)
{
  comm_ = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  map_ = Teuchos::rcp(new MapT(numGlobalUnknowns_,/*index base=*/0,comm_));
  numLocalUnknowns_ = static_cast<LO> (map_->getNodeNumElements()); 
  A_ = Teuchos::rcp(new MatrixT(map_,/*initialize to zero=*/true));
  b_ = Teuchos::rcp(new VectorT(map_,/*initialize to zero=*/true));
}

GO LinearSystem::mapLocalToGlobal(LO i)
{
  Teuchos::ArrayView<const GO> myGlobalElements =
    map_->getNodeElementList();
  return myGlobalElements[i];
}

void LinearSystem::insertToVector(double v, GO i)
{
  const ST value = static_cast<ST> (v);
  b_->replaceGlobalValue(i,value);
}

void LinearSystem::insertToMatrix(double v, GO i, GO j)
{
  const ST value = static_cast<ST> (v);
  A_->insertGlobalValues(i,
      Teuchos::tuple<GO>(j),
      Teuchos::tuple<ST>(value));
}

void LinearSystem::sumToVector(double v, GO i)
{
  const ST value = static_cast<ST> (v);
  b_->sumIntoGlobalValue(i,value);
}

void LinearSystem::sumToMatrix(double v, GO i, GO j)
{
  const ST value = static_cast<ST> (v);
  A_->sumIntoGlobalValues(i,
      Teuchos::tuple<GO>(j),
      Teuchos::tuple<ST>(value));
}

void LinearSystem::completeMatrixFill()
{
  A_->fillComplete();
}

void LinearSystem::describeVector()
{
  std::cout << b_->description() << std::endl;
}

void LinearSystem::describeMatrix()
{
  std::cout << A_->description() << std::endl;
}

void LinearSystem::solve()
{
}

}
