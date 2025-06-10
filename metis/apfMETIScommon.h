/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_METIS_COMMON_H
#define APF_METIS_COMMON_H

#include <vector>

namespace apf {

typedef NumberingOf<long> GlobalNumbering;
class Mesh;
class Migration;

namespace metis {

GlobalNumbering* makeNumbering(Mesh* m, const char* name, long start = 0);

void getOwnedAdjacencies(
  apf::GlobalNumbering* gn, std::vector<idx_t>& xadj,
  std::vector<idx_t>& adjncy, long gn_offset, bool remoteEdges
);

bool runMETIS(
  idx_t metis_nvtxs, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy,
  idx_t nparts, double imbalance, std::vector<idx_t>& part
);

apf::Migration* makePlan(
  GlobalNumbering* gn, const std::vector<idx_t>& owned_part, long gn_offset
);

} // namespace metis

} // namespace apf

#endif // APF_METIS_COMMON_H
