#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <pcu_util.h>

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
  // Localize n_owned_elms for Gatherv on xadj.
  // I reuse xadj_cts, xadj_displs below.
  std::vector<idx_t> part(mesh_->getPCU()->Self() == 0 ? metis_nvtxs : 0);
  idx_t nparts = multiple;
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
    lion_eprint(1, "ERROR: splitting failed: %s\n", metis_err);
    fail("metis splitting failed");
  }
  apf::Migration *plan = makePlan(gn, part, 0);
  apf::destroyGlobalNumbering(gn);
  auto t1 = pcu::Time();
  if (mesh_->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: split in %f seconds\n", t1 - t0);
  return plan;
}

} // namespace metis

} // namespace apf
