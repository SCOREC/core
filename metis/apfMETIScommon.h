/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <vector>

namespace apf {

typedef NumberingOf<long> GlobalNumbering;
class Mesh;
class Migration;

namespace metis {

void getOwnedAdjacencies(
  apf::GlobalNumbering* gn, std::vector<idx_t>& xadj,
  std::vector<idx_t>& adjncy, long gn_offset
);

apf::Migration* makePlan(
  GlobalNumbering* gn, const std::vector<idx_t>& owned_part, long gn_offset
);

} // namespace metis

} // namespace apf
