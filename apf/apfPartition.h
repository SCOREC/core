/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_PARTITION_H
#define APF_PARTITION_H

#include "apfMesh2.h"

namespace apf {

class Splitter
{
  public:
    virtual ~Splitter() {}
    virtual Migration* split(MeshTag* weights, double tolerance, int multiple) = 0;
};

class Balancer
{
  public:
    virtual ~Balancer() {}
    virtual void balance(MeshTag* weights, double tolerance) = 0;
};

struct Remap
{
  virtual int operator()(int n) = 0;
};

struct Divide : public Remap
{
  Divide(int n):by(n) {}
  int by;
  int operator()(int n) {return n / by;}
};

struct Multiply : public Remap
{
  Multiply(int n):by(n) {}
  int by;
  int operator()(int n) {return n * by;}
};

struct Modulo : public Remap
{
  Modulo(int n):by(n) {}
  int by;
  int operator()(int n) {return n % by;}
};

struct Unmodulo : public Remap
{
  Unmodulo(int original, int factor)
  {
    offset = (original / factor) * factor;
  }
  int offset;
  int operator()(int n) {return n + offset;}
};

struct Round : public Remap
{
  Round(int n):factor(n) {}
  int factor;
  int operator()(int n) {return (n / factor) * factor;}
};

void remapPartition(apf::Mesh2* m, Remap& remap);

}

#endif
