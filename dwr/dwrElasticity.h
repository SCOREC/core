/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_ELASTICITY_H
#define DWR_ELASTICITY_H

#include <apf.h>
#include <apfArray.h>
#include <apfDynamicArray.h>
#include <apfDynamicMatrix.h>
#include <Sacado.hpp>

namespace dwr {

typedef Sacado::Fad::DFad<double> AD;
typedef apf::Array<AD,3> AD_Vector3;
typedef apf::Array<AD_Vector3,3> AD_Matrix3x3;

class LinElastInt : public apf::Integrator
{
  public:
    LinElastInt(int o, apf::Field* u);
    ~LinElastInt();
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
