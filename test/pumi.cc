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
  pumi_printsys();

  double begin_mem = pumi_getmem(), begin_time=pumi_gettime();

  // read input args - in-model-file in-mesh-file out-mesh-file num-in-part
  getConfig(argc,argv);

  // load model
  gmi_model* g = pumi_geom_create(modelFile);
 
  // load mesh per process group
  int num_proc_group=1;
  if (num_in_part>1 && pumi_size()!=num_in_part)
    num_proc_group = pumi_size()/num_in_part;
  assert(pumi_size()%num_in_part==0);
  pMesh m=pumi_mesh_create(g, meshFile, num_in_part, num_proc_group);

  // write mesh in .smb
  pumi_mesh_write(m,outFile);
  // write mesh in .vtk
  pumi_mesh_write(m,"output", "vtk");

  // print mesh info
  int mesh_dim=pumi_mesh_getdim(m);
  if (!pumi_rank()) std::cout<<"[pumi_test]: mesh dim="<<mesh_dim<<"\n";

  // loop with mesh vertex
  int remote_count=0;
  pMeshEnt e;

  // loop over vertices
  pMeshIter it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    assert(pumi_ment_getdim(e)==0);
    assert(pumi_ment_getnumadj(e, mesh_dim+1)==0);
    if (!pumi_ment_isonbdry(e)) continue; // skip internal entity
    // if entity is on part boundary, count remote copies    
    pCopies copies;
    pumi_ment_getallrmt(e,copies);
    // loop over remote copies and increase the counter
    APF_ITERATE(pCopies,copies,rit)
      ++remote_count;
    // check #remotes
    assert (pumi_ment_getnumrmt(e)==copies.size());
    // check the entity is not ghost or ghosted
    assert(!pumi_ment_isghost(e) && !pumi_ment_isghosted(e));
  }
  m->end(it);

  if (!pumi_rank()) std::cout<<"[pumi_test]: remote_count="<<remote_count<<"\n";
  // print elapsed time and increased heap memory
  pumi_printtimemem("elapsed time and increased heap memory:", pumi_gettime()-begin_time, pumi_getmem()-begin_mem);

  // clean-up
  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
}

