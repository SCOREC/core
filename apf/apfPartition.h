/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_PARTITION_H
#define APF_PARTITION_H

#include "apfMesh.h"

namespace apf {

class Splitter
{
  public:
    virtual ~Splitter() {}
    virtual void split(MeshTag* weights, double tolerance, int multiple) = 0;
};

class Balancer
{
  public:
    virtual ~Balancer() {}
    virtual void balance(MeshTag* weights, double tolerance) = 0;
};

}

#endif
