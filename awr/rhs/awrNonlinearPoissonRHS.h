/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRNONLINEARPOISSONRHS_H
#define AWRNONLINEARPOISSONRHS_H

#include "awrRHS.h"

namespace awr {

class NonlinearPoissonRHS : public RHS
{
  public:
    NonlinearPoissonRHS(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~NonlinearPoissonRHS() {};
    virtual void
    evaluateElementRHS(apf::MeshEntity* element,
                       apf::DynamicMatrix& k);
  private:
    NonlinearPoissonRHS(const NonlinearPoissonRHS&);
    NonlinearPoissonRHS& operator=(const NonlinearPoissonRHS&);
};

}

#endif
