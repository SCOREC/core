#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>
#include <cassert>
#include <cstdlib>

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int num_in_part = 1;

void getConfig(int argc, char** argv)
{
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <num_in_mesh_part>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  num_in_part = atoi(argv[4]);
  assert(num_in_part <= PCU_Comm_Peers());
}

#include <pumi.h>
#include <iostream>
#include <cstdlib>
#include <assert.h>
#include <mpi.h>

int main_(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  if (argc<3)
  {
    if (!pumi_rank()) std::cout<<"[pumi_test] ./test_pumi geom_file.dmg mesh_file.smb num_in_part\n";
    pumi_finalize();
    return 1;
  }
  modelFile = argv[1];
  meshFile = argv[2];

  pGeom g = pumi_geom_create(modelFile);
  int num_in_part=1;
  if (argc>3 && atoi(argv[3])>1)
  {
    assert(atoi(argv[3]) == pumi_size());
    num_in_part == pumi_size();
  }
  pMesh m=pumi_mesh_create(g, meshFile, num_in_part);
  
  pumi_mesh_write(m, "output.smb", "mds");
  pumi_mesh_write(m, "output", "vtk");
  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
  return 1;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pumi_info();
  
  getConfig(argc,argv);

  gmi_model* g = pumi_geom_create(modelFile);
  pMesh m=pumi_mesh_create(g, meshFile, num_in_part);
  pumi_mesh_write(m,outFile);
  pumi_mesh_write(m,"output", "vtk");

  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
}

