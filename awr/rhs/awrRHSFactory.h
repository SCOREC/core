/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRRHSFACTORY_H
#define AWRRHSFACTORY_H

#include "awrRHS.h"
#include <Teuchos_RCP.hpp>

namespace awr {

class RHSFactory
{
  public:
    RHSFactory(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~RHSFactory() {};
    virtual Teuchos::RCP<RHS> create();
  protected:
    apf::Mesh* mesh_;
    Teuchos::ParameterList params_;
  private:
    RHSFactory(const RHSFactory&);
    RHSFactory& operator=(const RHSFactory&);
};

}

#endif
