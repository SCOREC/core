#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>
#include <cassert>
#include <cstdlib>
#include <pumi.h>

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int num_in_part = 0;

void getConfig(int argc, char** argv)
{
  if ( argc < 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <num_part_in_mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  if (argc>=4)
    num_in_part = atoi(argv[4]);
}

Ghosting* getGhostingPlan(pMesh m)
{
  int mesh_dim=m->getDimension();
  Ghosting* plan = new Ghosting(m, m->getDimension());
  if (!PCU_Comm_Self())
  {
  apf::MeshIterator* it = m->begin(mesh_dim);
  pMeshEnt e;
  int count=0;
  while ((e = m->iterate(it)))
  {
    int pid=(count+1)%PCU_Comm_Peers();
    plan->send(e, pid);
    if (count==10) break;
  }
  }
  return plan;
}

#include <iostream>
#include <cstdlib>
#include <mpi.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pumi_printsys();

#if 0
  int i, processid = getpid();
  if (!PCU_Comm_Self())
  {
    std::cout<<"Proc "<<PCU_Comm_Self()<<">> pid "<<processid<<" Enter any digit...\n";
    std::cin>>i;
  }
  else
    std::cout<<"Proc "<<PCU_Comm_Self()<<">> pid "<<processid<<" Waiting...\n";
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  double begin_mem = pumi_getmem(), begin_time=pumi_gettime();

  // read input args - in-model-file in-mesh-file out-mesh-file num-in-part
  getConfig(argc,argv);

  // load model
  gmi_model* g = pumi_geom_load(modelFile);
 
  // load mesh per process group
  assert(pumi_size()%num_in_part==0);
  if (!pumi_rank()) std::cout<<"[test_pumi] num_in_part="<<num_in_part<<"\n";
  pMesh m=NULL;
  if (num_in_part==1)
    m = pumi_mesh_loadserial(g, meshFile);
  else
    m = pumi_mesh_load(g, meshFile, num_in_part);

  pumi_mesh_print(m);
  pMeshEnt e;
  // let's do distribution
  // partitioning: sending an element to a single part
  // distribution: sending an element to multiple parts. Element may have remote copies.
  {
    Distribution* plan = new Distribution(m);

    int dim=pumi_mesh_getdim(m), count=0, pid;
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      pid=pumi_ment_getlocalid(e)%PCU_Comm_Peers();
      plan->send(e, pid);
      if (pid-1>=0) plan->send(e, pid-1);
      if (pid+1<PCU_Comm_Peers()) plan->send(e, pid+1);
      if (count==5) break;
      ++count;
    }
    m->end(it);
    plan->print(); // print distribution plan 
    pumi_mesh_distribute(m, plan);
  }  

  pumi_mesh_print(m);
  sleep(.5);
  if (!pumi_rank()) std::cout<<"\n";

  if (!pumi_rank()) std::cout<<"[test_pumi] writing mesh\n";

  // write mesh in .smb
  pumi_mesh_write(m,outFile);
  // write mesh in .vtk
  pumi_mesh_write(m,"output", "vtk");

  // print mesh info
  int mesh_dim=pumi_mesh_getdim(m);

  // loop with elements

  if (!pumi_rank()) std::cout<<"[test_pumi] checking various api's\n";
  pMeshIter mit = m->begin(mesh_dim);
  while ((e = m->iterate(mit)))
  {
    assert(pumi_ment_getdim(e)==mesh_dim);
    assert(pumi_ment_getnumadj(e, mesh_dim+1)==0);
    if (!pumi_ment_isonbdry(e)) continue; // skip internal entity
    // if entity is on part boundary, count remote copies    
    Copies copies;
    pumi_ment_getallrmt(e,copies);
    // loop over remote copies and increase the counter
    // check #remotes
    assert (pumi_ment_getnumrmt(e)==copies.size() && copies.size()>0);
    // check the entity is not ghost or ghosted
    assert(!pumi_ment_isghost(e) && !pumi_ment_isghosted(e));
  }
  m->end(mit);

  // re-load partitioned mesh
  {
    pumi_mesh_delete(m);
    g = pumi_geom_load(modelFile);
    m = pumi_mesh_load(g, meshFile, num_in_part); 
    pumi_mesh_print(m);
  }

  // let's do ghosting
  int num_org_vtx = pumi_mesh_getnument(m, 0);

  Ghosting* ghosting_plan = getGhostingPlan(m);
  ghosting_plan->print();
  //pumi_ghost_create(m, ghosting_plan);

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

  //pumi_mesh_verify(m);
 
  // print elapsed time and increased heap memory
  pumi_printtimemem("[test_pumi] elapsed time and increased heap memory:", pumi_gettime()-begin_time, pumi_getmem()-begin_mem);

  // clean-up
  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
}

