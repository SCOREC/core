/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWR_NONLINEARPOISSON_H
#define AWR_NONLINEARPOISSON_H

namespace awr {

class Problem;
class NonlinearPoissonIntegrator;

class NonlinearPoissonProblem : public Problem
{
  public:
    NonlinearPoissonProblem(ParameterList& p, apf::Mesh* m);
    ~NonlinearPoissonProblem();
  private:
    NonlinearPoissonIntegrator* integrator_;
    void validateProblemList();
    void setPrimalField();
    void createIntegrator();
    void processKe(apf::MeshEntity* e, apf::DynamicMatrix& Ke);
};

}

#endif
