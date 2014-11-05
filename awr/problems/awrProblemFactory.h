/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRPROBLEMFACTORY_H
#define AWRPROBLEMFACTORY_H

#include "awrProblem.h"

namespace awr {

class ProblemFactory
{
  public:
    ProblemFactory(apf::Mesh* m,
        const Teuchos::RCP<Teuchos::ParameterList>& p);
    virtual ~ProblemFactory() {};
    virtual Teuchos::RCP<Problem> create();
  protected:
    apf::Mesh* mesh_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
  private:
    ProblemFactory(const ProblemFactory&);
    ProblemFactory& operator=(const ProblemFactory&);
};

}

#endif
