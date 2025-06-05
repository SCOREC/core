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

#include <metis.h>

#include "apfMETISbalancer.h"
#include "apfMETIScommon.h"

namespace apf {

namespace metis {

static void gatherGraph(
  pcu::PCU& PCU, MPI_Comm comm,
  const std::vector<idx_t>& owned_xadj, const std::vector<idx_t>& owned_adjncy,
  std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy,
  std::vector<int>& vtx_cts, std::vector<int>& vtx_displs
) {
  auto t0 = pcu::Time();
  // Localize n_owned_elms for Gatherv on xadj.
  vtx_cts.resize(PCU.Self() == 0 ? PCU.Peers() : 0);
  int n_owned_elm = owned_xadj.size() - 1;
  MPI_Gather(
    &n_owned_elm, 1, MPI_INT,
    vtx_cts.data(), 1, MPI_INT,
    0, comm
  );
  if (PCU.Self() == 0) vtx_cts.back() += 1; // Add space for final xadj entry.
  vtx_displs.resize(PCU.Self() == 0 ? PCU.Peers() + 1 : 0);
  if (PCU.Self() == 0) {
    std::partial_sum(
      vtx_cts.begin(), vtx_cts.end(), vtx_displs.begin() + 1
    );
  }
  // Increment owned_xadj.
  int owned_adj_ct = owned_xadj.back();
  int xadj_offset = PCU.Exscan(owned_adj_ct);
  std::vector<idx_t> send_xadj(owned_xadj.size());
  for (size_t i = 0; i < owned_xadj.size(); ++i)
    send_xadj[i] = owned_xadj[i] + xadj_offset;
  // Gather xadj.
  xadj.resize(PCU.Self() == 0 ? vtx_displs.back() : 0);
  int xadj_sendct = send_xadj.size() - 1;
  // Final entry is equal to first of next rank, so only send on last rank.
  if (PCU.Self() == PCU.Peers() - 1) ++xadj_sendct;
  /* FIXME: use Pack/Unpack instead of direct MPI_Gatherv.
  PCU.Begin();
  PCU.Pack(0, xadj_sendct);
  PCU.Pack(0,
    owned_xadj.data(), xadj_sendct * sizeof(*owned_xadj.data())
  );
  PCU.Send();
  while (m->getPCU()->Receive()) {
    // Unpack...
  }
  */
  MPI_Gatherv(
    send_xadj.data(), xadj_sendct, MPI_INT,
    xadj.data(), vtx_cts.data(), vtx_displs.data(), MPI_INT,
    0, comm
  );
  // Undo "create the final xadj slot" from above.
  if (PCU.Self() == 0) vtx_cts.back() -= 1;
  // Gather adjncy.
  int owned_adjncy_ct = owned_adjncy.size();
  std::vector<int> adjncy_cts(PCU.Self() == 0 ? PCU.Peers() : 0);
  MPI_Gather(
    &owned_adjncy_ct, 1, MPI_INT,
    adjncy_cts.data(), 1, MPI_INT,
    0, comm
  );
  std::vector<int> adjncy_displs(PCU.Self() == 0 ? PCU.Peers() + 1 : 0);
  if (PCU.Self() == 0) {
    std::partial_sum(
      adjncy_cts.begin(), adjncy_cts.end(), adjncy_displs.begin() + 1
    );
  }
  adjncy.resize(PCU.Self() == 0 ? xadj.back() : 0);
  MPI_Gatherv(
    owned_adjncy.data(), owned_adjncy.size(), MPI_INT,
    adjncy.data(), adjncy_cts.data(), adjncy_displs.data(), MPI_INT,
    0, comm
  );
  auto t1 = pcu::Time();
  if (PCU.Self() == 0)
    lion_oprint(1, "METIS: localized graph in %f seconds\n", t1 - t0);
}

static void scatterPart(
  pcu::PCU& PCU, MPI_Comm comm,
  const std::vector<idx_t>& part,
  const std::vector<int>& vtx_cts, const std::vector<int>& vtx_displs,
  std::vector<idx_t>& owned_part, int n_owned
) {
  auto t0 = pcu::Time();
  owned_part.resize(n_owned);
  MPI_Scatterv(
    part.data(), vtx_cts.data(), vtx_displs.data(), MPI_INT,
    owned_part.data(), n_owned, MPI_INT, 0, comm
  );
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
  if (weights != nullptr) {
    if (mesh_->getPCU()->Self() == 0)
      lion_oprint(1, "METIS: ignoring weights\n");
  }
  // FIXME PCU_DEBUG_ASSERT(sizeof(idx_t) >= sizeof(mds_id_type));
  int elm_dim = mesh_->getDimension();
  PCU_ALWAYS_ASSERT(elm_dim == 3); // FIXME: update code to allow 2d
  auto t0 = pcu::Time();
  MPI_Comm comm;
  mesh_->getPCU()->DupComm(&comm);
  int n_owned_elm = apf::countOwned(mesh_, elm_dim);
  // Create global element numbering.
  auto numbering = apf::numberOwnedDimension(
    mesh_, "apfMETISbalancer_nb", elm_dim
  );
  // Don't use apf::makeGlobal because numbering may not be rank order.
  auto gn = apf::createGlobalNumbering(
    mesh_, "apfMETISbalancer_gnb", apf::getConstant(elm_dim)
  );
  int gn_offset = mesh_->getPCU()->Exscan(n_owned_elm);
  apf::MeshIterator *it = mesh_->begin(elm_dim);
  for (apf::MeshEntity *e; (e = mesh_->iterate(it));) {
    if (mesh_->isOwned(e)) {
      apf::number(gn, e, 0, gn_offset + apf::getNumber(numbering, e, 0, 0));
    }
  }
  mesh_->end(it);
  apf::synchronize(gn);
  std::vector<idx_t> owned_xadj, owned_adjncy;
  getOwnedAdjacencies(gn, owned_xadj, owned_adjncy, gn_offset, true);
  std::vector<idx_t> xadj, adjncy;
  std::vector<int> vtx_cts, vtx_displs;
  gatherGraph(
    *mesh_->getPCU(), comm, owned_xadj, owned_adjncy, xadj, adjncy,
    vtx_cts, vtx_displs
  );
  int metis_nvtxs = mesh_->getPCU()->Self() == 0 ? vtx_displs.back() - 1 : 0;
  std::vector<idx_t> part(mesh_->getPCU()->Self() == 0 ? metis_nvtxs : 0);
  if (mesh_->getPCU()->Self() == 0) {
    auto t0_part = pcu::Time();
    idx_t nparts = mesh_->getPCU()->Peers();
    std::vector<real_t> imb(nparts, tolerance);
    idx_t objval, ncon = 1;
    int r = METIS_PartGraphKway(
      &metis_nvtxs, &ncon, xadj.data(), adjncy.data(), // Graph
      NULL, NULL, NULL, // No sizing
      &nparts,
      NULL, imb.data(), NULL, &objval, part.data()
    );
    if (r != METIS_OK) {
      const char *metis_err = "";
      if (r == METIS_ERROR_INPUT) metis_err = "METIS: input error";
      else if (r == METIS_ERROR_MEMORY) metis_err = "METIS: memory error";
      else metis_err = "METIS: error";
      lion_eprint(1, "ERROR: balancing failed: %s\n", metis_err);
      return;
    }
    remapPart(nparts, part, vtx_cts);
    auto t1_part = pcu::Time();
    lion_oprint(1, "METIS: partitioned in %f seconds\n", t1_part - t0_part);
  }
  std::vector<idx_t> owned_part;
  scatterPart(
    *mesh_->getPCU(), comm, part, vtx_cts, vtx_displs, owned_part, n_owned_elm
  );
  apf::Migration *plan = makePlan(gn, owned_part, gn_offset);
  auto t0migrate = pcu::Time();
  mesh_->migrate(plan);
  auto t1migrate = pcu::Time();
  if (mesh_->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: migrated in %f seconds\n", t1migrate - t0migrate);
  plan = nullptr;
  apf::destroyGlobalNumbering(gn);
  apf::destroyNumbering(numbering);
  MPI_Comm_free(&comm);
  auto t1 = pcu::Time();
  if (mesh_->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: balanced in %f seconds\n", t1 - t0);
}

} // namespace metis

} // namespace apf

