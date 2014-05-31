/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ZOLTAN_H
#define APF_ZOLTAN_H

namespace apf {

enum {
  RCB,
  RIB,
  HYPERGRAPH,
  PARMETIS,
  GRAPH
};

enum {
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

Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach,
    bool debug = true, bool sync = true);
Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach,
    bool debug = true);

}

#endif
