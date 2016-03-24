#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <cassert>
#include <cstdlib>
#include <iostream>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int partitionFactor = 1;

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

void switchToOriginals()
{
  int self = PCU_Comm_Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}

void getConfig(int argc, char** argv)
{
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
  assert(partitionFactor <= PCU_Comm_Peers());
}

}

apf::Migration* get_model_plan(gmi_model* g, apf::Mesh* m)
{
  apf::Migration* plan = new apf::Migration(m);
  apf::MeshEntity* e;
  int num_gface = g->n[2];
  assert(gmi_find(g, 2, num_gface));

  int gface_id;
  gmi_ent* gface;
  apf::MeshIterator* it = m->begin(2);
  while ((e = m->iterate(it)))
  {
    gface = (gmi_ent*)(m->toModel(e));
    gface_id = gmi_tag(g, gface);
    plan->send(e, gface_id-1);
  }
  m->end(it);
  return plan;
}

#include <cstring>
#include "omega_h.h"
#include <apfOmega_h.h>
#include "apfMesh2.h"
#include "apf.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  getConfig(argc,argv);
  bool isOriginal = ((PCU_Comm_Self() % partitionFactor) == 0);
  gmi_model* g = 0;
  g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  switchToOriginals();
  if (isOriginal) {
    m = apf::loadMdsMesh(g, meshFile);
    plan = get_model_plan(g, m);
  }
  switchToAll();
  m = repeatMdsMesh(m, g, plan, partitionFactor);
  m->writeNative(outFile);

  // write to vtk
  apf::writeVtkFiles("partitioned", m);

  // ghosting
  apf::Mesh2* gm = apf::makeEmptyMdsMesh(g, 2, false);

  double t0 = PCU_Time();
  osh_t osh_mesh = osh::fromAPF(m);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    std::cout << "to omega_h " << t1 - t0 << " seconds\n";
  osh_ghost(osh_mesh, 3);
  double t2 = PCU_Time();
  if (!PCU_Comm_Self())
    std::cout << "3 layer ghosting " << t2 - t1 << " seconds\n";
  osh::toAPF(osh_mesh, gm);
  double t3 = PCU_Time();
  if (!PCU_Comm_Self())
    std::cout << "back to APF " << t3 - t2 << " seconds\n";

  // attach ghost flag

  apf::MeshTag* own_tag = gm->findTag("owner");
  assert (own_tag);
  int own_partid;

  apf::Field* ghost_f = apf::createStepField(gm, "ghost_field", apf::SCALAR);
  apf::Field* own_f = apf::createStepField(gm, "own_field", apf::SCALAR);
  apf::MeshIterator* it = gm->begin(2);
  apf::MeshEntity* e;
  double ghost_value;
  int ghost_flag=1, non_ghost_flag=0;
  while ((e = gm->iterate(it)))
  {
    gm->getIntTag(e, own_tag, &own_partid); //gm->getOwner(e) does not work
    if (own_partid!=PCU_Comm_Self())
      ghost_value=ghost_flag;
    else
      ghost_value=non_ghost_flag;
    setScalar(ghost_f, e, 0, ghost_value);
    setScalar(own_f, e, 0, own_partid);
  }
  gm->end(it);

  apf::writeVtkFiles("ghosted", gm);

  destroyField(ghost_f);
  destroyField(own_f);
  osh_free(osh_mesh);

  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

