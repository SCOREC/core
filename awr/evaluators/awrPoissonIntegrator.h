/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apf.h>
#include <apfShape.h>
#include <apfDynamicMatrix.h>

namespace awr {

class PoissonIntegrator : public apf::Integrator
{
  public:
    PoissonIntegrator(apf::Field* f, int o);
    void inElement(apf::MeshElement* me);
    void outElement(apf::MeshElement* e) {};
    void atPoint(apf::Vector3 const &p, double w, double dv);
    void parallelReduce() {};
  private:
    apf::Field* f_;
    apf::FieldShape* shape_;
    apf::MeshElement* me_;
    apf::Element* e_;
    int numComponents_;
    int numNodes_; 
    int numDims_;
    int numLocalEqs_;
    apf::DynamicMatrix Ke_;
};

}
