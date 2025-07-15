/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMETIS.h"
#include "apfMETISbalancer.h"
#include "apfMETISsplitter.h"

namespace apf {

Splitter* makeMETISsplitter(Mesh* mesh) {
  return new metis::MetisSplitter(mesh);
}

Balancer* makeMETISbalancer(Mesh* mesh) {
  return new metis::MetisBalancer(mesh);
}

bool hasMETIS() { return true; }

} // namespace apf

