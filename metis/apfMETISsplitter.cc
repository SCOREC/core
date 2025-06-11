#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include "apfMETIS.h"

#include <metis.h>

#include "apfMETISsplitter.h"
#include "apfMETIScommon.h"

namespace apf {

namespace metis {

Migration* MetisSplitter::split(
  MeshTag* weights, double tolerance, int multiple
) {
  PCU_ALWAYS_ASSERT(tolerance > 1.0);
  PCU_ALWAYS_ASSERT(multiple > 1);
  if (mesh_->getPCU()->Peers() > APF_METIS_MAXRANKS) {
    fail(
      "METIS called with > " STRINGIFY(APF_METIS_MAXRANKS)
      " procs, which is unsupported due to memory requirements\n"
    );
  }
  auto t0 = pcu::Time();
  if (weights != nullptr) {
    lion_oprint(1, "METIS: weights are not supported\n");
  }
  int elm_dim = mesh_->getDimension();
  PCU_ALWAYS_ASSERT(elm_dim == 3); // FIXME: update code to allow 2d
  int metis_nvtxs = apf::countOwned(mesh_, elm_dim);
  auto gn = makeNumbering(mesh_, "apfMETISsplitter_gnb");
  std::vector<idx_t> xadj, adjncy;
  getOwnedAdjacencies(gn, xadj, adjncy, 0, false);
  std::vector<idx_t> part;
  bool r = runMETIS(metis_nvtxs, xadj, adjncy, multiple, tolerance, part);
  if (!r) fail("METIS splitting failed");
  apf::Migration *plan = makePlan(gn, part, 0);
  apf::destroyGlobalNumbering(gn);
  auto t1 = pcu::Time();
  if (mesh_->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: split in %f seconds\n", t1 - t0);
  // Spread parts over PCU::Peers() * factor parts.
  if (mesh_->getPCU()->Self() != 0) {
    for (int i = 0; i < plan->count(); ++i) {
      apf::MeshEntity *e = plan->get(i);
      int dest = plan->sending(e);
      dest += mesh_->getPCU()->Self() * multiple;
      plan->send(e, dest);
    }
  }
  return plan;
}

} // namespace metis

} // namespace apf
