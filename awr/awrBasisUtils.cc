/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrBasisUtils.h"

namespace awr {

BasisUtils::BasisUtils(apf::Mesh* m) :
  mesh_(m)
{
}

int BasisUtils::getNumDims()
{
}

int BasisUtils::getNumNodes(apf::Field* f)
{
}

int BasisUtils::getNumQP()
{
}

void BasisUtils::getBF(apf::Element* elem, NodeQPScalar& bf)
{
}

void BasisUtils::getWBF(apf::Element* elem, NodeQPScalar& w_bf)
{
}

void BasisUtils::getGradBF(apf::Element* elem, NodeQPVector& grad_bf)
{
}

void BasisUtils::getWGradBF(apf::Element* elem, NodeQPVector& w_grad_bf)
{
}

}
