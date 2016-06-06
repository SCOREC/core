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
int num_in_part = 0;
int num_proc_group = 1;

void getConfig(int argc, char** argv)
{
  if ( argc < 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <num_in_mesh_part>  <num_proc_group>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  if (argc>4)
    num_in_part = atoi(argv[4]);
  else 
    num_in_part=PCU_Comm_Peers();
  if (argc>5) 
    num_proc_group = atoi(argv[5]);
  assert(num_in_part <= PCU_Comm_Peers());
  if (argc==5 && num_in_part!=1)
    num_proc_group = PCU_Comm_Peers()/num_in_part;

}

#include <pumi.h>
#include <iostream>
#include <cstdlib>
#include <mpi.h>

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
  pumi_mesh_print(m);
  sleep(.5);
  if (!pumi_rank()) std::cout<<"\n";

  // write mesh in .smb
  pumi_mesh_write(m,outFile);
  // write mesh in .vtk
  pumi_mesh_write(m,"output", "vtk");

  // print mesh info
  int mesh_dim=pumi_mesh_getdim(m);

  // loop with mesh vertex
  int remote_count=0;
  pMeshEnt e;

  // loop over vertices
  pMeshIter mit = m->begin(0);
  while ((e = m->iterate(mit)))
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
  m->end(mit);
  int num_org_vtx = pumi_mesh_getnument(m, 0);

  m = pumi_ghost_create(0, mesh_dim, 2, 1);
  if (!pumi_rank()) std::cout<<"creating ghost layers\n";
  int num_ghost_vtx=0;
  mit = m->begin(0);
  while ((e = m->iterate(mit)))
  {
    if (pumi_ment_isghost(e))
    {
      ++num_ghost_vtx;
     assert(pumi_ment_getownpid(e)!=pumi_rank());
    }
  }   
  m->end(mit);
  assert(num_ghost_vtx+num_org_vtx==pumi_mesh_getnument(m,0));

  // print mesh info
  pumi_mesh_print(m);
  if (!pumi_rank()) std::cout<<"\n";

  // FIXME: deleting ghost layers is temporarily unavailable
  //pumi_ghost_delete(m);
 
  // print elapsed time and increased heap memory
  pumi_printtimemem("elapsed time and increased heap memory:", pumi_gettime()-begin_time, pumi_getmem()-begin_mem);

  // clean-up
  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
}

