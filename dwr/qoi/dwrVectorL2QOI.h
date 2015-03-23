/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_VECTORL2QOI_H
#define DWR_VECTORL2QOI_H

#include <apf.h>
#include <apfDynamicVector.h>

namespace dwr {

class VectorL2QOI : public apf::Integrator
{
  public:
    VectorL2QOI(int o, apf::Field* u);
    ~VectorL2QOI();
    void inElement(apf::MeshElement* me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    apf::DynamicVector Fe;
  private:
    int numDims_;
    int numNodes_;
    int numDofs_;
    apf::Element* e_;
    apf::Field* u_;
};

}

#endif
