/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRLHSFACTORY_H
#define AWRLHSFACTORY_H

#include "awrLHS.h"
#include <Teuchos_RCP.hpp>

namespace awr {

class LHSFactory
{
  public:
    LHSFactory(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~LHSFactory() {};
    virtual Teuchos::RCP<LHS> create();
  protected:
    apf::Mesh* mesh_;
    Teuchos::ParameterList params_;
  private:
    LHSFactory(const LHSFactory&);
    LHSFactory& operator=(const LHSFactory&);
};

}

#endif
