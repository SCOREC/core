#include <apfMesh.h>
#include <apfShape.h>
#include <lionPrint.h>
#include "maBalance.h"
#include "maAdapt.h"
#include <parma.h>
#include <apfZoltan.h>
#if defined(APW_LGMETIS)
#include <numeric>
#include "maDBG.h"
#endif

#define MAX_ZOLTAN_GRAPH_RANKS 16*1024

namespace ma {

static double clamp(double x, double max, double min)
{
  if (x > max) return max;
  if (x < min) return min;
  return x;
}

static double getSizeWeight(Adapt* a, Entity* e, int type)
{
/* until we have a size field that matches the
   anisotropy of the layer itself, we have to
   hack around to prevent prisms from getting
   very small weights due to their very small thickness.
   the current hack will be to measure their
   triangular bases instead.
   pyramids will get very small weights, but
   they are supposed to be fairly sparse so
   that should not ruin things. */
  if (type == apf::Mesh::PRISM) {
    Entity* f[5];
    a->mesh->getDownward(e, 2, f);
    return a->sizeField->getWeight(f[0]);
  }
  return a->sizeField->getWeight(e);
}

static double clampForIterations(Adapt* a, double weight)
{
  Mesh* m = a->mesh;
  int dimension = m->getDimension();
  double max = pow(2.0, dimension*(a->refinesLeft));
/* coarsening performance is more empirical: 3x decrease in tet
   count when uniformly refining a 58k element cube, 4x decrease
   on some 2D meshes which were more structured. */
  double min = pow(4.0, -(a->coarsensLeft));
  return clamp(weight, max, min);
}

static double clampForLayerPermissions(Adapt* a, int type, double weight)
{
  if (apf::isSimplex(type))
    return weight;
  if ( ! a->input->shouldRefineLayer)
    weight = std::max(1.0, weight);
  if ( ! a->input->shouldCoarsenLayer)
    weight = std::min(1.0, weight);
  return weight;
}

static double accountForTets(Adapt* a, int type, double weight)
{
  if (a->input->shouldTurnLayerToTets) {
    if (type == apf::Mesh::PRISM)
      return weight * 3;
    if (type == apf::Mesh::PYRAMID)
      return weight * 2;
  }
  return weight;
}

double getElementWeight(Adapt* a, Entity* e)
{
  int type = a->mesh->getType(e);
  double weight = getSizeWeight(a, e, type);
  weight = clampForIterations(a, weight);
  weight = clampForLayerPermissions(a, type, weight);
  return accountForTets(a, type, weight);
}

Tag* getElementWeights(Adapt* a)
{
  Mesh* m = a->mesh;
  Tag* weights = m->createDoubleTag("ma_weight",1);
  Entity* e;
  Iterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
  {
    double weight = getElementWeight(a,e);
    m->setDoubleTag(e,weights,&weight);
  }
  m->end(it);
  return weights;
}

static void runBalancer(Adapt* a, apf::Balancer* b)
{
  Mesh* m = a->mesh;
  Input* in = a->input;
  Tag* weights = getElementWeights(a);
  b->balance(weights,in->maximumImbalance);
  delete b;
  removeTagFromDimension(m,weights,m->getDimension());
  m->destroyTag(weights);
}

void runZoltan(Adapt* a, int method=apf::GRAPH)
{
  runBalancer(a, apf::makeZoltanBalancer(
        a->mesh, method, apf::REPARTITION,
        /* debug = */ false));
}

void runParma(Adapt* a)
{
  runBalancer(a, Parma_MakeElmBalancer(a->mesh));
}

void printEntityImbalance(Mesh* m)
{
  double imbalance[4];
  Parma_GetEntImbalance(m,&imbalance);
  double p = (imbalance[m->getDimension()]-1)*100;
  print(m->getPCU(), "element imbalance %.0f%% of average", p);
}

double estimateWeightedImbalance(Adapt* a)
{
  Tag* w = getElementWeights(a);
  double imb[4];
  Parma_GetWeightedEntImbalance(a->mesh, w, &imb);
  removeTagFromDimension(a->mesh, w, a->mesh->getDimension());
  a->mesh->destroyTag(w);
  return imb[a->mesh->getDimension()];
}

#ifdef APW_LGMETIS
#include <metis.h>

void runLocalizedGraphMetis(Adapt* a) {
  auto t0 = pcu::Time();
  MPI_Comm comm;
  a->mesh->getPCU()->DupComm(&comm);
  int elm_dim = a->mesh->getDimension();
  // FIXME: confirm sizeof(idx_t) >= sizeof(mds_id_type)
  int nelm = apf::countOwned(a->mesh, elm_dim);
  std::vector<int> nelms(a->mesh->getPCU()->Peers());
  a->mesh->getPCU()->Allgather(&nelm, nelms.data(), 1);
  std::vector<int> nelm_sums(nelms.size() + 1);
  nelm_sums[0] = 0;
  std::partial_sum(nelms.begin(), nelms.end(), nelm_sums.begin() + 1);
  idx_t nvtxs = nelm_sums.back(), ncon = 1;
  // Create global element numbering.
  auto numbering = apf::numberOwnedDimension(
    a->mesh, "maBalance__runLGM_nb", elm_dim
  );
  // Don't use apf::makeGlobal because it will recalculate nelm_sums.
  auto gn = apf::createGlobalNumbering(
    a->mesh, "maBalance_runLGM_gnb", apf::getConstant(elm_dim)
  );
  int gn_offset = nelm_sums[a->mesh->getPCU()->Self()];
  apf::MeshIterator *it = a->mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = a->mesh->iterate(it));) {
    if (a->mesh->isOwned(e)) {
      apf::number(gn, e, 0, gn_offset + apf::getNumber(numbering, e, 0, 0));
    }
  }
  a->mesh->end(it);
  apf::synchronize(gn);
  auto opposites = apf::tagOpposites(gn, "maBalance__run_LGM_opp");
  // Build local adj_cts.
  std::vector<idx_t> owned_adj_cts(nelm);
  it = a->mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = a->mesh->iterate(it));) {
    int e_num = apf::getNumber(numbering, e, 0, 0);
    apf::Downward graph_edges;
    int nd = a->mesh->getDownward(e, elm_dim - 1, graph_edges);
    int nd_own = 0;
    for (int i = 0; i < nd; ++i) {
      // If mesh face (3d)/edge (2d) is shared or has 2 upward adjacencies,
      // then add a graph edge slot.
      if (a->mesh->isShared(graph_edges[i])) ++nd_own;
      else if (a->mesh->countUpward(graph_edges[i]) == 2) ++nd_own;
    }
    owned_adj_cts[e_num] = nd_own;
  }
  a->mesh->end(it);
  // Build local xadj FIXME: inefficient
  std::vector<idx_t> owned_xadj(1 + owned_adj_cts.size());
  std::partial_sum(
    owned_adj_cts.begin(), owned_adj_cts.end(), owned_xadj.begin() + 1
  );
  // Build local adjncy.
  std::vector<idx_t> owned_adjncy(owned_xadj.back());
  it = a->mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = a->mesh->iterate(it));) {
    int e_num = apf::getNumber(numbering, e, 0, 0);
    int e_xadj = owned_xadj[e_num];
    apf::Downward graph_edges;
    int nd = a->mesh->getDownward(e, elm_dim - 1, graph_edges);
    int adj_i = 0;
    for (int j = 0; j < nd; ++j) {
      if (a->mesh->isShared(graph_edges[j])) {
        long opp_num;
        a->mesh->getLongTag(graph_edges[j], opposites, &opp_num);
        owned_adjncy[e_xadj + adj_i] = opp_num;
        ++adj_i;
      } else {
        PCU_DEBUG_ASSERT(a->mesh->countUpward(graph_edges[j]) <= 2);
        for (int k = 0; k < a->mesh->countUpward(graph_edges[j]); ++k) {
          apf::MeshEntity *up = a->mesh->getUpward(graph_edges[j], k);
          if (up != e) {
            owned_adjncy[e_xadj + adj_i] = apf::getNumber(gn, up, 0);
            ++adj_i;
          }
        }
      }
    }
  }
  a->mesh->end(it);
  // Build xadj.
  idx_t owned_adjncy_size = owned_adjncy.size();
  std::vector<idx_t> owned_adjncy_sizes(a->mesh->getPCU()->Peers());
  a->mesh->getPCU()->Allgather(
    &owned_adjncy_size, owned_adjncy_sizes.data(), 1
  );
  std::vector<idx_t> adj_cts(nvtxs);
  MPI_Gatherv(
    owned_adj_cts.data(), owned_adj_cts.size(), MPI_INT,
    adj_cts.data(), nelms.data(), nelm_sums.data(), MPI_INT,
    0, comm
  );
  std::vector<idx_t> xadj(nvtxs + 1);
  std::partial_sum(adj_cts.begin(), adj_cts.end(), xadj.begin() + 1);
  std::vector<idx_t> adjncy(xadj.back());
  MPI_Gatherv(
    owned_adjncy.data(), owned_adjncy.size(), MPI_INT,
    adjncy.data(), owned_adjncy_sizes.data(), xadj.data(), MPI_INT,
    0, comm
  );
  auto t_localize = pcu::Time();
  if (a->mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: localized graph in %f seconds\n", t_localize - t0);
  std::vector<idx_t> part(nvtxs);
  if (a->mesh->getPCU()->Self() == 0) {
#ifdef APW_LGMETIS_SER
    idx_t nparts = 4; // a->mesh->getPCU()->Peers();
#else
    idx_t nparts = a->mesh->getPCU()->Peers();
#endif
    std::vector<real_t> imb(nparts, a->input->maximumImbalance);
    idx_t objval;
    int r = METIS_PartGraphKway(
      &nvtxs, &ncon, xadj.data(), adjncy.data(), // Graph
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
  }
  auto t_part = pcu::Time();
  if (a->mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: partitioned in %f seconds\n", t_part - t_localize);
  std::vector<idx_t> owned_part(nelm);
  MPI_Scatterv(
    part.data(), nelms.data(), nelm_sums.data(), MPI_INT,
    owned_part.data(), nelm, MPI_INT, 0, comm
  );
  auto t_scatter = pcu::Time();
  if (a->mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: scattered in %f seconds\n", t_scatter - t_part);
#if defined(APW_LGMETIS_VIZ)
  apf::Field* dst_part = apf::createStepField(a->mesh, "dst_part", apf::SCALAR);
#endif
  apf::Migration *plan = new apf::Migration(a->mesh);
  it = a->mesh->begin(elm_dim);
  for (apf::MeshEntity *e; (e = a->mesh->iterate(it));) {
    int dest = owned_part[apf::getNumber(gn, e, 0) - gn_offset];
    if (dest != a->mesh->getPCU()->Self()) plan->send(e, dest);
#if defined(APW_LGMETIS_VIZ)
    apf::setScalar(dst_part, e, 0, dest);
#endif
  }
  auto t_plan = pcu::Time();
  if (a->mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: planned in %f seconds\n", t_plan - t_scatter);
#if defined(APW_LGMETIS_VIZ)
  static int num_vtk = 0;
  ma_dbg::dumpMeshWithQualities(a, num_vtk, "migrate");
  ++num_vtk;
  apf::destroyField(dst_part);
#endif
  auto t0migrate = pcu::Time();
#ifdef APW_LGMETIS_SER
  delete plan;
#else
  a->mesh->migrate(plan);
#endif
  auto t1migrate = pcu::Time();
  if (a->mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: migrated in %f seconds\n", t1migrate - t0migrate);
  plan = nullptr;
  apf::removeTagFromDimension(a->mesh, opposites, elm_dim - 1);
  a->mesh->destroyTag(opposites);
  apf::destroyGlobalNumbering(gn);
  apf::destroyNumbering(numbering);
  MPI_Comm_free(&comm);
  auto t1 = pcu::Time();
  if (a->mesh->getPCU()->Self() == 0)
    lion_oprint(1, "METIS: balanced in %f seconds\n", t1 - t0);
#if defined(APW_LGMETIS_VIZ)
  it = a->mesh->begin(elm_dim - 1);
  int shared = 0;
  for (apf::MeshEntity* e; (e = a->mesh->iterate(it));) {
    if (a->mesh->isShared(e)) ++shared;
  }
  a->mesh->end(it);
  lion_oprint(1, "METIS: rank %d has %d shared faces\n",
    a->mesh->getPCU()->Self(), shared
  );
#endif
}
#endif

void preBalance(Adapt* a)
{
#ifndef APW_LGMETIS_SER
  if (a->mesh->getPCU()->Peers()==1)
    return;
#endif
#ifdef APW_LGMETIS
  runLocalizedGraphMetis(a);
  return;
#endif
  Input* in = a->input;
  // First take care of user overrides. That is, if any of the three options
  // is true, apply that balancer and return.
  if (in->shouldRunPreZoltan) {
    runZoltan(a);
    return;
  }
  if (in->shouldRunPreZoltanRib) {
    runZoltan(a,apf::RIB);
    return;
  }
  if (in->shouldRunPreParma) {
    runParma(a);
    return;
  }

  // Then, take care of the case where all the options are set to false.
  // That is, if the default values have not changed by the user. In
  // this case, we apply the best possible balancer, if weighted imbalance
  // is bigger than in->maximumImbalance
  if ((!in->shouldRunPreZoltan) &&
      (!in->shouldRunPreZoltanRib) &&
      (!in->shouldRunPreParma) &&
      (estimateWeightedImbalance(a) > in->maximumImbalance)) {
#ifdef PUMI_HAS_ZOLTAN
    // The parmetis multi-level graph partitioner memory usage grows
    // significantly with process count beyond 16K processes
    if (a->mesh->getPCU()->Peers() < MAX_ZOLTAN_GRAPH_RANKS) {
      runZoltan(a);
      return;
    }
    else {
      runZoltan(a, apf::RIB);
      return;
    }
#else
    runParma(a);
    return;
#endif
  }
}

void midBalance(Adapt* a)
{
  if (a->mesh->getPCU()->Peers()==1)
    return;
#ifdef APW_LGMETIS
  runLocalizedGraphMetis(a);
  return;
#else
  Input* in = a->input;
  // First take care of user overrides. That is, if any of the three options
  // is true, apply that balancer and return.
  if (in->shouldRunMidZoltan) {
    runZoltan(a);
    return;
  }
  if (in->shouldRunMidParma) {
    runParma(a);
    return;
  }
  // Then, take care of the case where all the options are set to false.
  // That is, if the default values have not changed by the user. In
  // this case, we apply the best possible balancer, if weighted imbalance
  // is bigger than in->maximumImbalance
  if ((!in->shouldRunMidZoltan) &&
      (!in->shouldRunMidParma) &&
      (estimateWeightedImbalance(a) > in->maximumImbalance)) {
#ifdef PUMI_HAS_ZOLTAN
    // The parmetis multi-level graph partitioner memory usage grows
    // significantly with process count beyond 16K processes
    if (a->mesh->getPCU()->Peers() < MAX_ZOLTAN_GRAPH_RANKS) {
      runZoltan(a);
      return;
    }
    else {
      runZoltan(a, apf::RIB);
      return;
    }
#else
    runParma(a);
    return;
#endif
  }
#endif
}

void postBalance(Adapt* a)
{
  if (a->mesh->getPCU()->Peers()==1)
    return;
#ifdef APW_LGMETIS
  runLocalizedGraphMetis(a);
  return;
#endif
  Input* in = a->input;
  // First take care of user overrides. That is, if any of the three options
  // is true, apply that balancer and return.
  if (in->shouldRunPostZoltan) {
    runZoltan(a);
    printEntityImbalance(a->mesh);
    return;
  }
  if (in->shouldRunPostZoltanRib) {
    runZoltan(a,apf::RIB);
    printEntityImbalance(a->mesh);
    return;
  }
  if (in->shouldRunPostParma) {
    runParma(a);
    printEntityImbalance(a->mesh);
    return;
  }
  // Then, take care of the case where all the options are set to false.
  // That is, if the default values have not changed by the user. In
  // this case, we apply the best possible balancer, if weighted imbalance
  // is bigger than in->maximumImbalance
  if ((!in->shouldRunPostZoltan) &&
      (!in->shouldRunPostZoltanRib) &&
      (!in->shouldRunPostParma) &&
      (estimateWeightedImbalance(a) > in->maximumImbalance)) {
#ifdef PUMI_HAS_ZOLTAN
    // The parmetis multi-level graph partitioner memory usage grows
    // significantly with process count beyond 16K processes
    if (a->mesh->getPCU()->Peers() < MAX_ZOLTAN_GRAPH_RANKS) {
      runZoltan(a);
      printEntityImbalance(a->mesh);
      return;
    }
    else {
      runZoltan(a, apf::RIB);
      printEntityImbalance(a->mesh);
      return;
    }
#else
    runParma(a);
    printEntityImbalance(a->mesh);
    return;
#endif
  }
}

}
