/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef AWRBASISUTILS_H
#define AWRBASISUTILS_H

#include "apf.h"
#include "apfDynamicArray.h"

namespace awr {

typedef
apf::DynamicArray<
apf::DynamicArray<double> > NodeQPScalar;

typedef
apf::DynamicArray<
apf::DynamicArray<
apf::DynamicArray<double> > > NodeQPVector;


class BasisUtils
{
  public:
    BasisUtils(apf::Mesh* m);
    int getNumDims();
    int getNumNodes(apf::Field* f);
    int getNumQP();
    void getBF(apf::Element* elem, NodeQPScalar& bf);
    void getWBF(apf::Element* elem, NodeQPScalar& w_bf);
    void getGradBF(apf::Element* elem, NodeQPVector& grad_bf);
    void getWGradBF(apf::Element* elem, NodeQPVector& w_grad_bf);
  private:
    apf::Mesh* mesh_;
};

}

#endif
