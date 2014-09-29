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
    BasisUtils(apf::Mesh* m, apf::Field* f, apf::MeshEntity* e, int o);
    ~BasisUtils();
    int getNumDims();
    int getNumQP();
    int getNumNodes();
    void getBF(NodeQPScalar& bf);
    void getWBF(NodeQPScalar& w_bf);
    void getGradBF(NodeQPVector& grad_bf);
    void getWGradBF(NodeQPVector& w_grad_bf);
  private:
    apf::Mesh* mesh_;
    apf::Field* sol_;
    apf::MeshEntity* elem_;
    apf::MeshElement* mesh_elem_;
    apf::Element* field_elem_;
    int order_;
};

}

#endif
