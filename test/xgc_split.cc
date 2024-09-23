#include "pumi.h"

#include <apf.h>
#include <cstring>
#include <mpi.h>
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>
#include "apfMDS.h"
#include "apfShape.h"

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int serial=0;

void getConfig(int argc, char** argv, pcu::PCU* PCUObj)
{
  if (argc < 4) {
    if (!PCUObj->Self())
      printf("Usage: %s <model> <mesh> <outMesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  if (argc > 4) 
    serial=atoi(argv[4]);
}

Migration* get_xgc_plan(pGeom g, pMesh m)
{
  int dim = pumi_mesh_getDim(m);
  Migration* plan = new Migration(m);
  if (!pumi_rank()) return plan;

  pMeshEnt e;
  int num_gface = pumi_geom_getNumEnt(g, dim);
  PCU_ALWAYS_ASSERT(num_gface==pumi_size());
  int gface_id;
  int dest_pid;
  pMeshIter it = m->begin(2); // face
  while ((e = m->iterate(it))) 
  { 
    pGeomEnt gface = pumi_ment_getGeomClas(e); // get the classification
    gface_id = pumi_gent_getID(gface); // get the geom face id
    dest_pid = gface_id-1;
    plan->send(e, dest_pid);    
  }
  m->end(it);
  return plan;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  pumi_load_pcu(&PCUObj);
  getConfig(argc,argv,&PCUObj);

  pGeom g = pumi_geom_load(modelFile);
  pMesh m;
  if (serial) 
  {
    m = pumi_mesh_loadSerial(g, meshFile);
    // split a serial mesh based on model ID
    Migration* plan = get_xgc_plan(g, m);
    pumi_mesh_migrate(m, plan);
    pumi_mesh_write(m, outFile);
  }
  else 
    m = pumi_mesh_load(g, meshFile, pumi_size());

  // write to vtk
  char without_extension[256];
  snprintf(without_extension,strlen(argv[3])-3,"%s",argv[3]);

  char vtk_fname[512];
  snprintf(vtk_fname,512,"%s",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");

  // ghosting
  pumi_ghost_createLayer(m, 0, 2, 3, 0);
  snprintf(vtk_fname,512,"%s-ghosted",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");
  pumi_ghost_delete(m);

  pumi_ghost_createLayer(m, 0, 2, 3, 1);
  snprintf(vtk_fname,512,"%s-ghosted-copy",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");
  pumi_ghost_delete(m);

  pumi_mesh_delete(m);

  }
  MPI_Finalize();
}

