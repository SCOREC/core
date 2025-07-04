/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include "apfMETIS.h"

#include <metis.h>

#include "apfMETISbalancer.h"
#include "apfMETIScommon.h"

namespace apf {

namespace metis {

static void gatherGraph(
  pcu::PCU& PCU,
  std::vector<idx_t>& owned_xadj, const std::vector<idx_t>& owned_adjncy,
  std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy,
  std::vector<int>& vtx_cts
) {
  PCU_DEBUG_ASSERT(PCU.Peers() > 1);
  PCU_DEBUG_ASSERT(!owned_xadj.empty());
  PCU_DEBUG_ASSERT(owned_adjncy.size() == size_t(owned_xadj.back()));
  auto t0 = pcu::Time();
  int owned_vtx_ct = owned_xadj.size() - 1;
  int xadj_size = PCU.Add(owned_vtx_ct) + 1;
  xadj.resize(PCU.Self() == 0 ? xadj_size : 0);
  vtx_cts.resize(PCU.Self() == 0 ? PCU.Peers() : 0);
  // Increment owned_xadj.
  int xadj_offset = PCU.Exscan(owned_xadj.back());
  for (size_t i = 0; i < owned_xadj.size(); ++i) owned_xadj[i] += xadj_offset;
  int xadj_sendct = owned_vtx_ct;
  int xadj_displ = PCU.Exscan(xadj_sendct);
  // Final owned_xadj entry is equal to first entry on the next rank, so only
  // send it on the last rank:
  if (PCU.Self() == PCU.Peers() - 1) ++xadj_sendct;
  // Gather xadj.
  PCU.Begin();
  int pcu_self = PCU.Self();
  PCU.Pack(0, pcu_self);
  PCU.Pack(0, owned_vtx_ct);
  PCU.Pack(0, xadj_displ);
  PCU.Pack(0, owned_xadj.data(), xadj_sendct * sizeof(*owned_xadj.data()));
  PCU.Send();
  while (PCU.Receive()) {
    int rank, displ;
    PCU.Unpack(rank);
    PCU.Unpack(vtx_cts[rank]);
    PCU.Unpack(displ);
    int recvct = vtx_cts[rank];
    if (rank == PCU.Peers() - 1) ++recvct;
    PCU.Unpack(xadj.data() + displ, recvct * sizeof(*xadj.data()));
  }
  // Gather adjncy.
  adjncy.resize(PCU.Self() == 0 ? xadj.back() : 0);
  int adjncy_ct = owned_adjncy.size();
  int adjncy_displ = PCU.Exscan(adjncy_ct);
  PCU.Begin();
  PCU.Pack(0, adjncy_ct);
  PCU.Pack(0, adjncy_displ);
  PCU.Pack(0, owned_adjncy.data(), adjncy_ct * sizeof(*adjncy.data()));
  PCU.Send();
  while (PCU.Receive()) {
    int ct, displ;
    PCU.Unpack(ct);
    PCU.Unpack(displ);
    PCU.Unpack(adjncy.data() + displ, ct * sizeof(*adjncy.data()));
  }
  auto t1 = pcu::Time();
  if (PCU.Self() == 0)
    lion_oprint(1, "METIS: localized graph in %f seconds\n", t1 - t0);
}

static void scatterPart(
  pcu::PCU& PCU,
  const std::vector<idx_t>& part,
  const std::vector<int>& vtx_cts,
  std::vector<idx_t>& owned_part, int n_owned
) {
  PCU_DEBUG_ASSERT(PCU.Peers() > 1);
  PCU_DEBUG_ASSERT(PCU.Self() != 0 || vtx_cts.size() == size_t(PCU.Peers()));
  PCU_DEBUG_ASSERT(
    std::accumulate(vtx_cts.begin(), vtx_cts.end(), 0UL) == part.size()
  );
  auto t0 = pcu::Time();
  owned_part.resize(n_owned);
  PCU.Begin();
  if (PCU.Self() == 0) {
    int displ = 0;
    for (size_t i = 0; i < vtx_cts.size(); ++i) {
      PCU.Pack(i, part.data() + displ, vtx_cts[i] * sizeof(*part.data()));
      displ += vtx_cts[i];
    }
  }
  PCU.Send();
  while (PCU.Receive()) {
    PCU.Unpack(owned_part.data(), n_owned * sizeof(*owned_part.data()));
  }
  auto t1 = pcu::Time();
  if (PCU.Self() == 0)
    lion_oprint(1, "METIS: scattered in %f seconds\n", t1 - t0);
}

static void remapPart(int nparts, std::vector<idx_t>& part, const std::vector<int>& owned_cts) {
  std::map<int,int> dest_map; // map from metis_part to dst_part.
  // Loop through each part (maybe randomize order)
  for (int p = 0, i = 0; p < nparts; ++p) {
    // Count elements for each destination part.
    std::map<int,int> dest_counts;
    for (int j = 0; j < owned_cts[p]; ++i, ++j) {
      int dest = part[i];
      ++dest_counts[dest];
    }
    // Remap max not in dest_map to original part.
    while (!dest_counts.empty()) {
      using dest_counts_v = decltype(dest_counts)::value_type;
      auto max_it = std::max_element(dest_counts.cbegin(), dest_counts.cend(),
        [](const dest_counts_v& a, const dest_counts_v& b) {
          return a.second < b.second;
        }
      );
      int dest = max_it->first;
      if (dest_map.count(dest)) dest_counts.erase(max_it);
      else {
        // Update set mapped.
        dest_map[dest] = p;
        dest_counts.clear();
      }
    }
  }
  // Confirm remap is valid.
  bool valid_map = true;
  for (int p = 0; p < nparts && valid_map; ++p) {
    if (dest_map.count(p) != 1) valid_map = false;
  }
  if (!valid_map) {
    lion_oprint(1, "METIS: fixing up mapping\n");
    // Collect parts not mapped TO.
    std::set<int> unmapped;
    for (int p = 0; p < nparts; ++p) unmapped.insert(p);
    for (auto p : dest_map) unmapped.erase(p.second);
    // Iterate through parts not mapped FROM and assign mapping.
    for (int p = 0; p < nparts; ++p) {
      if (dest_map.count(p) != 1) {
        dest_map[p] = *unmapped.begin();
        unmapped.erase(unmapped.begin());
      }
    }
  }
  // Remap destination parts.
  for (int p = 0, i = 0; p < nparts; ++p) {
    for (int j = 0; j < owned_cts[p]; ++i, ++j) part[i] = dest_map[part[i]];
  }
}

void MetisBalancer::balance(MeshTag* weights, double tolerance) {
  PCU_ALWAYS_ASSERT(tolerance > 1.0);
  if (mesh_->getPCU()->Peers() == 1) return; // no work to be done.
  else if (mesh_->getPCU()->Peers() > APF_METIS_MAXRANKS) {
    fail(
      "METIS called with > " STRINGIFY(APF_METIS_MAXRANKS)
      " procs, which is unsupported due to memory requirements\n"
    );
  }
  if (weights != nullptr) {
    if (mesh_->getPCU()->Self() == 0)
      lion_oprint(1, "METIS: weights are not supported\n");
  }
  int elm_dim = mesh_->getDimension();
  PCU_ALWAYS_ASSERT(elm_dim == 3); // FIXME: update code to allow 2d
  auto t0 = pcu::Time();
  int n_owned_elm = apf::countOwned(mesh_, elm_dim);
  long gn_offset = mesh_->getPCU()->Exscan(long(n_owned_elm));
  auto gn = makeNumbering(mesh_, "apfMETISbalancer_gnb", gn_offset);
  std::vector<idx_t> owned_xadj, owned_adjncy;
  getOwnedAdjacencies(gn, owned_xadj, owned_adjncy, gn_offset, true);
  std::vector<idx_t> xadj, adjncy;
  std::vector<int> vtx_cts;
  gatherGraph(
    *mesh_->getPCU(), owned_xadj, owned_adjncy, xadj, adjncy, vtx_cts
  );
  int metis_nvtxs = std::accumulate(vtx_cts.begin(), vtx_cts.end(), 0);
  std::vector<idx_t> part;
  if (mesh_->getPCU()->Self() == 0) {
    idx_t nparts = mesh_->getPCU()->Peers();
    bool r = runMETIS(metis_nvtxs, xadj, adjncy, nparts, tolerance, part);
    if (!r) {
      lion_eprint(1, "ERROR: balancing failed. continuning...\n");
      return;
    }
    auto t0_remap = pcu::Time();
    remapPart(nparts, part, vtx_cts);
    auto t1_remap = pcu::Time();
    lion_oprint(1, "METIS: remapped in %f seconds\n", t1_remap - t0_remap);
  }
  std::vector<idx_t> owned_part;
  scatterPart(
    *mesh_->getPCU(), part, vtx_cts, owned_part, n_owned_elm
  );
  apf::Migration *plan = makePlan(gn, owned_part, gn_offset);
  auto t0migrate = pcu::Time();
  mesh_->migrate(plan);
  auto t1migrate = pcu::Time();
  if (mesh_->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: migrated in %f seconds\n", t1migrate - t0migrate);
  plan = nullptr;
  apf::destroyGlobalNumbering(gn);
  auto t1 = pcu::Time();
  if (mesh_->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: balanced in %f seconds\n", t1 - t0);
}

} // namespace metis

} // namespace apf

