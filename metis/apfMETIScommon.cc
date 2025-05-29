/*
 * Copyright (C) 2025 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <numeric>

#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfZoltan.h> // FIXME: tagOpposites
#include <lionPrint.h>
#include <pcu_util.h>

#include <metis.h>

namespace apf {

namespace metis {

void getOwnedAdjacencies(
  GlobalNumbering* gn, std::vector<idx_t>& xadj,
  std::vector<idx_t>& adjncy, long gn_offset, bool remoteEdges
) {
  Mesh* mesh = apf::getMesh(gn);
  int elm_dim = mesh->getDimension();
  int n_owned_elm = apf::countOwned(mesh, elm_dim);
  MeshTag* opp_tag = nullptr;
  if (remoteEdges) {
    opp_tag = apf::tagOpposites(gn, "maBalance__run_LGM_opp");
  }
  // Build local adj_cts.
  std::vector<idx_t> adj_cts(n_owned_elm);
  MeshIterator* it = mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = mesh->iterate(it));) {
    int e_num = apf::getNumber(gn, e, 0) - gn_offset;
    apf::Downward graph_edges;
    int nd = mesh->getDownward(e, elm_dim - 1, graph_edges);
    int nd_own = 0;
    for (int i = 0; i < nd; ++i) {
      // If mesh face (3d)/edge (2d) is shared or has 2 upward adjacencies,
      // then add a graph edge slot.
      if (mesh->isShared(graph_edges[i])) {
        if (remoteEdges) ++nd_own;
      } else if (mesh->countUpward(graph_edges[i]) == 2) ++nd_own;
    }
    PCU_DEBUG_ASSERT(e_num < adj_cts.size());
    adj_cts[e_num] = nd_own;
  }
  mesh->end(it);
  // Build local xadj.
  xadj.resize(n_owned_elm + 1);
  std::partial_sum(
    adj_cts.begin(), adj_cts.end(), xadj.begin() + 1
  );
  // Build local adjncy.
  adjncy.resize(xadj.back());
  it = mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = mesh->iterate(it));) {
    int local_num = apf::getNumber(gn, e, 0) - gn_offset;
    PCU_DEBUG_ASSERT(local_num < xadj.size());
    int e_xadj = xadj[local_num]; // FIXME: double check.
    apf::Downward graph_edges;
    int nd = mesh->getDownward(e, elm_dim - 1, graph_edges);
    int adj_i = 0;
    for (int j = 0; j < nd; ++j) {
      if (mesh->isShared(graph_edges[j])) {
        if (remoteEdges) {
          long opp_num;
          mesh->getLongTag(graph_edges[j], opp_tag, &opp_num);
          PCU_DEBUG_ASSERT(e_xadj + adj_i < adjncy.size());
          adjncy[e_xadj + adj_i] = opp_num;
          ++adj_i;
        }
      } else {
        // FIXME: assert may not be true in non-manifold 2d cases.
        //        the loop may still make sense.
        PCU_DEBUG_ASSERT(mesh->countUpward(graph_edges[j]) <= 2);
        for (int k = 0; k < mesh->countUpward(graph_edges[j]); ++k) {
          apf::MeshEntity *up = mesh->getUpward(graph_edges[j], k);
          if (up != e) {
            PCU_DEBUG_ASSERT(e_xadj + adj_i < adjncy.size());
            adjncy[e_xadj + adj_i] = apf::getNumber(gn, up, 0);
            ++adj_i;
          }
        }
      }
    }
  }
  mesh->end(it);
  if (remoteEdges) {
    apf::removeTagFromDimension(mesh, opp_tag, elm_dim - 1);
    mesh->destroyTag(opp_tag);
  }
}

apf::Migration* makePlan(
  GlobalNumbering* gn, const std::vector<idx_t>& owned_part, long gn_offset
) {
  Mesh* mesh = getMesh(gn);
  int elm_dim = mesh->getDimension();
  auto t0 = pcu::Time();
  apf::Migration *plan = new apf::Migration(mesh);
  MeshIterator* it = mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = mesh->iterate(it));) {
    long local_num = apf::getNumber(gn, e, 0) - gn_offset;
    PCU_DEBUG_ASSERT(local_num < owned_part.size());
    int dest = owned_part[local_num];
    if (dest != mesh->getPCU()->Self()) plan->send(e, dest);
  }
  auto t1 = pcu::Time();
  if (mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: planned in %f seconds\n", t1 - t0);
  return plan;
}

} // namespace metis

} // namespace apf
