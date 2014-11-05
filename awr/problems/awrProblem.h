/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRPROBLEM_H
#define AWRPROBLEM_H

#include "../linear_system/awrLinearSystem.h"
#include <apf.h>

namespace awr {

class Problem
{
  public:
    Problem(apf::Mesh* m, const Teuchos::RCP<Teuchos::ParameterList>& p);
    virtual ~Problem() {};
  protected:
    apf::Mesh* mesh_;
    int numEquations_;
    GO numGlobalUnknowns_;
    Teuchos::RCP<LinearSystem> linearSystem_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
  private:
    void validateSublists();
    Problem(const Problem&);
    Problem& operator=(const Problem&);
};

}

#endif
