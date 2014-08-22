/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ZOLTAN_H
#define APF_ZOLTAN_H

/** \file apfZoltan.h
  \brief Zoltan partitioning for apf::Mesh objects */

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
  PARTITION,
  REPARTITION,
  REFINE,
  PART_KWAY, 
  PART_GEOM,
  PART_GEOM_KWAY,
  ADAPT_REPART,
  REFINE_KWAY
};

class Mesh;
class Splitter;
class Balancer;

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

}

#endif
