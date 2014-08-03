#include <PCU.h>
#include "parma.h"
#include <mpi.h>
#include <apfMesh2.h>
#include <apfMDS.h> /* yea, implementation dependent for now. */

static bool amChosen(int factor)
{
  return (PCU_Comm_Self() % factor) == 0;
}

static int nearestChosen(int factor)
{
  int self = PCU_Comm_Self();
  return self - (self % factor);
}

static void migrateToChosen(apf::Mesh2* m, int factor)
{
  apf::Migration* plan = new apf::Migration(m);
  if ( ! amChosen(factor)) {
    int to = nearestChosen(factor);
    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    while ((e = m->iterate(it)))
      plan->send(e, to);
    m->end(it);
  }
  if (!PCU_Comm_Self())
    fprintf(stderr,"migrating to multiples of %d (ignore empty part warning)\n",
        factor);
  m->migrate(plan);
}

static MPI_Comm shrinkComm(int factor)
{
  int self = PCU_Comm_Self();
  int newSelf = self / factor;
  int group = self % factor;
  MPI_Comm newComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, newSelf, &newComm);
  return newComm;
}

static apf::Migration* finishChosen(apf::Mesh2* m, int factor)
{
  apf::Splitter* s = Parma_MakeRibSplitter(m);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  double const doesntMatter = 1.10;
  apf::Migration* plan = s->split(weights, doesntMatter, factor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  return plan;
}

static void expand(apf::Mesh2* m, int factor, apf::Migration* plan)
{
  if ( ! amChosen(factor))
    plan = new apf::Migration(m);
  m->migrate(plan);
}

namespace alb {

void shrinkPartition(apf::Mesh2* m, int factor, void (*runAfter)(apf::Mesh2* m))
{
  migrateToChosen(m, factor);
  MPI_Comm subComm = shrinkComm(factor);
  bool chosen = amChosen(factor);
  shrinkMdsPartition(m, factor);
  PCU_Switch_Comm(subComm);
  apf::Migration* plan = 0;
  if (chosen) {
    runAfter(m);
    plan = finishChosen(m, factor);
    expandMdsPartition(m, factor);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&subComm);
  expand(m, factor, plan);
}

}
