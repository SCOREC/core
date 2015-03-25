/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_ELASTICITYLHS_H
#define DWR_ELASTICITYLHS_H

#include <apf.h>
#include <apfDynamicArray.h>
#include <apfDynamicMatrix.h>
#include "dwrTypes.h"

namespace dwr {

class ElasticityLHS : public apf::Integrator
{
  public:
    ElasticityLHS(int o, apf::Field* u);
    ~ElasticityLHS();
    void inElement(apf::MeshElement* me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    void setElasticModulus(double E) { E_ = E; }
    void setPoissonsRatio(double nu) { nu_ = nu; }
    apf::DynamicMatrix Ke;
  private:
    double E_;
    double nu_;
    int numDims_;
    int numNodes_;
    int numDofs_;
    apf::Element* e_;
    apf::Field* primal_;
    apf::DynamicArray<AD_Vector3> u_;
    apf::DynamicArray<AD> residual_;
};

}

#endif
