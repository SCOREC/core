#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <cstring>
#include "pumi.h"
#include <PCU.h>
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
  if ( argc != 4) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = pumi_size();
  assert(partitionFactor <= PCU_Comm_Peers());
}

}

apf::Migration* get_xgc_plan(gmi_model* g, apf::Mesh* m, int num_peers)
{
  int dim = m->getDimension();
  apf::Migration* plan = new apf::Migration(m);
  apf::MeshEntity* e;
  int num_gface = g->n[2];
  gmi_ent* gface = gmi_find(g, 2, num_gface);
  assert(gface);
  int gface_id;
  int dest_pid;
  apf::MeshIterator* it = m->begin(2); // face
  while ((e = m->iterate(it))) 
  { 
    gface = (gmi_ent*)(m->toModel(e)); // get the classification
    gface_id = gmi_tag(g, gface); // get the geom face id
    dest_pid = gface_id-1;
    plan->send(e, dest_pid);    
  }
  m->end(it);
  return plan;
}

pMesh split_serial_xgc_mesh(pGeom model)
{
  bool isOriginal = ((PCU_Comm_Self() % partitionFactor) == 0);

  gmi_model* g = model->getGmi();

  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
//  apf::MeshTag* dest_pid_tag = m->createIntTag("dest_pid", 2);
  switchToOriginals();
  if (isOriginal) {
    m = pumi_mesh_load(model, meshFile, 1);
    plan = get_xgc_plan(g, m, pumi_size()); // determine the destination part id for elements
  }
  switchToAll();
  m = repeatMdsMesh(m, g, plan, partitionFactor);
  pumi_mesh_write(m,outFile);

  // write to vtk
  char without_extension[256];
  snprintf(without_extension,strlen(outFile)-3,"%s",outFile);
  char vtk_fname[32];
  sprintf(vtk_fname,"%s.vtk",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");
  m->verify();
  return m;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();

  getConfig(argc,argv);

  pGeom model = pumi_geom_load(modelFile);
  // FIXME: loading serial xgc mesh (/meshes/xgc/mesh.smb) then splitting and ghosting fails inside ghosting
  //  pMesh m = split_serial_xgc_mesh(model);
  // Loading a partitioned xgc mesh (/meshes/xgc/1-xpoint.smb) and ghosting works
  pMesh m = pumi_mesh_load(model, meshFile, pumi_size());

  // ghosting
  pumi_ghost_createLayer(m, 0, 2, 1, 1);

  // attach ghost flag
  apf::Field* ghost_f = apf::createStepField(m, "ghost_field", apf::SCALAR);
  apf::Field* own_f = apf::createStepField(m, "own_field", apf::SCALAR);
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  double ghost_value, own_partid;
  int ghost_flag=1, non_ghost_flag=0;
  while ((e = m->iterate(it))) 
  { 
    own_partid=pumi_ment_getOwnPID(e);
    if (pumi_ment_isGhost(e))
      ghost_value=ghost_flag;
    else
      ghost_value=non_ghost_flag;
    setScalar(ghost_f, e, 0, ghost_value);
    setScalar(own_f, e, 0, own_partid);
  }
  m->end(it);

  pumi_mesh_write(m, "ghosted", "vtk");

  pumi_mesh_delete(m);

  pumi_finalize();
  MPI_Finalize();
}

