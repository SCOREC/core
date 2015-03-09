/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "dwr.h"

namespace dwr {

ElasticityProblem* createElasticityProblem()
{
  return new ElasticityProblem();
}

void destroyElasticityProblem(ElasticityProblem* p)
{
  delete p;
}

}
