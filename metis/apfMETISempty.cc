/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apf.h>
#include "apfMETIS.h"

namespace apf {

Splitter* makeMETISsplitter(Mesh*) {
  fail("apf_metis compiled without METIS support!");
  return 0;
}

Balancer* makeMETISbalancer(Mesh*) {
  fail("apf_metis compiled without METIS support!");
  return 0;
}

bool hasMETIS() { return false; }

} // namespace apf
