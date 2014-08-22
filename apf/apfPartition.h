/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_PARTITION_H
#define APF_PARTITION_H

/** \file apfPartition.h
  \brief tools for changing mesh partitioning */

#include "apfMesh2.h"

namespace apf {

/** \brief Splits a mesh part into many
 \details usually the construction of a specific Splitter will
          ask for a mesh object, which is the part to be acted upon */
class Splitter
{
  public:
    virtual ~Splitter() {}
    /** \brief call the underlying split algorithm
      \param weights a tag of one double that should be attached to all the
                     elements. splitters will try to balance the total element
                     weight of the resulting parts.
                     Some splitters support a weights pointer of zero, in
                     which case elements are weighted equally
      \param tolerance a factor greater than one at such that if
                       max_weight < tolerance * average_weight
                       the splitter will consider the parts balanced
      \param multiple how many output parts to generate
      \returns a migration object mapping all elements in the local part
               to part ids from 0 to multiple - 1.
               some splitters have convenience options that will change
               the resulting part ids. */
    virtual Migration* split(MeshTag* weights, double tolerance,
                             int multiple) = 0;
};

/** \brief Load balance over all mesh parts
 \details usually the construction of a specific Balancer will
          ask for a mesh object, which is the mesh to be acted upon */
class Balancer
{
  public:
    virtual ~Balancer() {}
    /** \brief call collective load balancing
      \param weights a tag of one double that should be attached to all the
                     elements. balancers will try to balance the total element
                     weight of the existing parts.
                     Some splitters support a weights pointer of zero, in
                     which case elements are weighted equally
      \param tolerance a factor greater than one at such that if
                       max_weight < tolerance * average_weight
                       the balancer will consider the parts balanced
      \details the Balancer is responsible for running migrations
               to improve the load balance by its implementation-defined
               criteria. It is not allowed to modify anything besides
               partitioning */
    virtual void balance(MeshTag* weights, double tolerance) = 0;
};

/** \brief a map from old part ids to new part ids */
struct Remap
{
  virtual int operator()(int n) = 0;
};

/** \brief divide the part id */
struct Divide : public Remap
{
  Divide(int n):by(n) {}
  int operator()(int n) {return n / by;}
private:
  int by;
};

/** \brief multiply the part id */
struct Multiply : public Remap
{
  Multiply(int n):by(n) {}
  int operator()(int n) {return n * by;}
private:
  int by;
};

/** \brief return part id modulo n */
struct Modulo : public Remap
{
  Modulo(int n):by(n) {}
  int operator()(int n) {return n % by;}
private:
  int by;
};

/** \brief inverse of apf::Modulo */
struct Unmodulo : public Remap
{
  Unmodulo(int original, int factor)
  {
    offset = (original / factor) * factor;
  }
  int operator()(int n) {return n + offset;}
private:  
  int offset;
};

/** \brief map to nearest multiple of n */
struct Round : public Remap
{
  Round(int n):factor(n) {}
  int operator()(int n) {return (n / factor) * factor;}
private:
  int factor;
};

/** \brief remap all part ids in the mesh structure
  \details when using sub-group partitioning schemes
           or splitting meshes (see Parma_SplitPartition or apf::Splitter),
           it is useful to be able to update all partition model
           structures in a mesh to reflect a transition from one
           partitioning scheme to the next.

           this function applies the given map to all part ids in the
           remote copies, resident sets, and matching using the apf::Mesh2
           interface */
void remapPartition(apf::Mesh2* m, Remap& remap);

}

#endif
