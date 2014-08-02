#include <PCU.h>
#include "parma.h"
#include <mpi.h>
#include <apfMesh2.h>

static bool amChosen(int factor)
{
  return (PCU_Comm_Self() % factor) == 0;
}

static bool nearestChosen(int factor)
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
  m->migrate(plan);
}

static MPI_Comm shrinkComm(int factor)
{
  int self = PCU_Comm_Self();
  int group = self % factor;
  MPI_Comm newComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, self, &newComm);
  return newComm;
}

static apf::Migration* finishChosen(apf::Mesh2* m, int factor)
{
  apf::Splitter* s = Parma_MakeRibSplitter(m);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = s->split(weights, 1.10, factor);
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

void Parma_ShrinkPartition(apf::Mesh2* m, int factor,
    void (*runAfter)(apf::Mesh2* m))
{
  migrateToChosen(m, factor);
  MPI_Comm subComm = shrinkComm(factor);
  bool chosen = amChosen(factor);
  PCU_Switch_Comm(subComm);
  apf::Migration* plan = 0;
  if (chosen) {
    runAfter(m);
    plan = finishChosen(m, factor);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  PCU_Switch_Comm(MPI_COMM_WORLD);
  expand(m, factor, plan);
}
