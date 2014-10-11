/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRNONLINEARPOISSONLHS_H
#define AWRNONLINEARPOISSONLHS_H

#include "awrLHS.h"

namespace awr {

class NonlinearPoissonLHS : public LHS
{
  public:
    NonlinearPoissonLHS(apf::Mesh* m, const Teuchos::ParameterList& p);
    virtual ~NonlinearPoissonLHS() {};
    virtual void
    evaluateElementLHS(apf::MeshEntity* element,
                       apf::DynamicMatrix& k);
  private:
    NonlinearPoissonLHS(const NonlinearPoissonLHS&);
    NonlinearPoissonLHS& operator=(const NonlinearPoissonLHS&);
};

}

#endif
