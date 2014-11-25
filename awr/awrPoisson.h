/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_POISSON_H
#define AWR_POISSON_H

namespace awr {

class Problem;
class PoissonIntegrator;

class PoissonProblem : public Problem
{
  public:
    PoissonProblem(ParameterList& p, apf::Mesh* m);
    ~PoissonProblem();
  private:
    PoissonIntegrator* integrator_;
    void validateProblemList();
    void setPrimalField();
    void createIntegrator();
    void processKe(apf::MeshEntity* e, apf::DynamicMatrix& Ke);
};

}

#endif
