/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRPOISSONPROBLEM_H
#define AWRPOISSONPROBLEM_H

#include "awrProblem.h"

namespace awr {

class PoissonProblem : public Problem
{
  public:
    PoissonProblem(apf::Mesh* m,
        const Teuchos::RCP<Teuchos::ParameterList>& p);
    virtual ~PoissonProblem() {};
  private:
    int integrationOrder_;
    apf::Field* primalSolution_;
    void validateParameters();
    PoissonProblem(const PoissonProblem&);
    PoissonProblem& operator=(const PoissonProblem&);
};

}

#endif
