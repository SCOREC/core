/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ZOLTAN_H
#define APF_ZOLTAN_H

/** \page zoltan APF-Zoltan

  Zoltan is a unification of several scientific data partitioning tools
  which are especially suited to parallel mesh partitioning.
  It is developed at Sandia National Labs, and its home page is
  http://www.cs.sandia.gov/zoltan/

  This package implements an interface between APF and Zoltan,
  Converting data from an apf::Mesh to Zoltan data structures,
  running one of the many algorithms available in Zoltan, and
  using the result to migrate the mesh.
  By using a mesh abstraction together with a partitioner abstraction,
  we create a very general and useful component.

  The interface is in apfZoltan.h
  */

/** \file apfZoltan.h
  \brief Zoltan partitioning for apf::Mesh objects */

#include <apfNumbering.h>

namespace apf {

/** \brief Zoltan partitioning method */
enum ZoltanMethod {
  /** \brief Recursive Coordinate Bisection */
  RCB,
  /** \brief Recursive Inertial Bisection */
  RIB,
  /** \brief Hyper-graph partitioning */
  HYPERGRAPH,
  /** \brief Use ParMetis */
  PARMETIS,
  /** \brief General graph partitionig */
  GRAPH
};

/** \brief Zoltan partitioning approach */
enum ZoltanApproach {
  /** \brief (Hyper)Graph - does not consider the initial distribution */
  PARTITION,
  /** \brief (Hyper)Graph - considers the initial distribution */
  REPARTITION,
  /** \brief (HYPER)Graph - targets partitions needing only small changes */
  REFINE,
  /** \brief Graph - multilevel */
  PART_KWAY,
  /** \brief Graph - space filling curves */
  PART_GEOM,
  /** \brief Graph - hybrid method combining PART_KWAY and PART_GEOM */
  PART_GEOM_KWAY,
  /** \brief Graph - targets graphs generated from adaptively refined meshes */
  ADAPT_REPART,
  /** \brief Graph - targets partitions needing only small changes*/
  REFINE_KWAY
};

class Mesh;
class Splitter;
class Balancer;
class MeshTag;

/** \brief Make a Zoltan Splitter object
  \details the resulting splitter will apply Zoltan
  to the local mesh part to break it into several new parts.
  \param method select from apf::ZoltanMethod
  \param approach select from apf::ZoltanApproach
  \param debug print the full Zoltan configuration when splitting
  \param sync all parts are splitting by the same factor,
              multiply the part ids in the resulting apf::Migration
              accordingly */
Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach,
    bool debug = true, bool sync = true);

/** \brief Make a Zoltan Splitter object
  \details the resulting splitter will apply Zoltan
  to the global mesh part to break it into several new parts.
  \param method select from apf::ZoltanMethod
  \param approach select from apf::ZoltanApproach
  \param debug print the full Zoltan configuration when splitting
  */
Splitter* makeZoltanGlobalSplitter(Mesh* mesh, int method, int approach,
    bool debug = true);

/** \brief Make a Zoltan Balancer object
  \details this Balancer will apply Zoltan to the global mesh
  to improve load balance.

  Also note that this Balancer will create a Zoltan edge
  between two elements that share matched faces.
  \param method select from apf::ZoltanMethod
  \param approach select from apf::ZoltanApproach
  \param debug print the full Zoltan configuration */
Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach,
    bool debug = true);

/** \brief Tag global ids of opposite elements to boundary faces
  \details this function creates a LONG tag of one value
  and attaches to all partition boundary faces the global
  id of the element on the other side.
  \param gn global element numbering
  \param name the name of the resulting tag */
MeshTag* tagOpposites(GlobalNumbering* gn, const char* name);

/** \brief Get an element-to-element connectivity array
  \details this function assumes the mesh has one element type.
  the resulting array is created with new int[nelements * nsides].
  nsides is the number of faces of an element.
  entry [i * nsides + j] is the global id of the j'th adjacent
  element to local element i, which can be -1 for a geometric
  boundary. */
int* getElementToElement(apf::Mesh* m);

}

#endif
