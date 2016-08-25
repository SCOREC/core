#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfZoltan.h>
#include <cassert>
#include <cstdlib>
#include <pumi.h>
#include <unistd.h>
#include <iostream>
#include <mpi.h>

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int num_in_part = 0;
int do_distr=0;

void getConfig(int argc, char** argv)
{
  if ( argc < 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <num_part_in_mesh> <do_distribution(0/1)>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  if (argc>=4)
    num_in_part = atoi(argv[4]);
if (argc>=5)
    do_distr = atoi(argv[5]);
}

Ghosting* getGhostingPlan(pMesh m)
{
  int mesh_dim=m->getDimension();
  Ghosting* plan = new Ghosting(m, m->getDimension());
  {
    apf::MeshIterator* it = m->begin(mesh_dim);
    pMeshEnt e;
    int count=0;
    while ((e = m->iterate(it)))
    {
      int pid = (pumi_ment_getglobalid(e)+1)%pumi_size();
      plan->send(e, pid);
      ++count; 
     if (count==m->count(mesh_dim)/5) break;
    }
    m->end(it);
  }
  return plan;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pumi_printsys();

#if 1
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



  // read input args - in-model-file in-mesh-file out-mesh-file num-in-part
  getConfig(argc,argv);

  // load model
  gmi_model* g = pumi_geom_load(modelFile);
 
  // load mesh per process group
  assert(pumi_size()%num_in_part==0);
  if (!pumi_rank()) std::cout<<"[test_pumi] num_in_part="<<num_in_part<<"\n";

  double begin_mem = pumi_getmem(), begin_time=pumi_gettime();

  pMesh m=NULL;
  if (do_distr)
    m = pumi_mesh_loadserial(g, meshFile);
  else
  {
    m = pumi_mesh_load(g, meshFile, num_in_part); // static partitioning if num_in_part=1
    if (num_in_part==1) 
      pumi_mesh_write(m,"mesh.smb");
  }

  pMeshEnt e;

  // distribution: sending an element to multiple parts. Element may have remote copies.
  if (do_distr)
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

    pumi_mesh_distribute(m, plan);
    if (!pumi_rank()) std::cout<<"\n[test_pumi] writing mesh in vtk\n";

    // write mesh in .smb
    pumi_mesh_write(m,outFile);  
  }  

  pumi_mesh_write(m,"output", "vtk");

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

  if (1)// re-load partitioned mesh
  {
    pumi_mesh_delete(m);
    g = pumi_geom_load(modelFile);
    if (num_in_part==1 && pumi_size()>1)
      m =  pumi_mesh_load(g, "mesh.smb", pumi_size()); 
    else
      m = pumi_mesh_load(g, meshFile, num_in_part); 
  }
  if (!pumi_rank()) std::cout<<"\n[test_pumi] delete and reload mesh\n";
  pumi_mesh_print(m);

  // let's do ghosting
  int num_org_vtx = pumi_mesh_getnument(m, 0);
  int* org_mcount=new int[4];
  for (int i=0; i<4; ++i)
  {
    org_mcount[i] = m->count(i);
    std::cout<<"("<<pumi_rank()<<") INFO dim "<<i<<": org ent count "<<org_mcount[i]<<"\n";
  }

  Ghosting* ghosting_plan = getGhostingPlan(m);

  pumi_ghost_create(m, ghosting_plan);

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

  pumi_ghost_delete(m);

  for (int i=0; i<4; ++i)
    assert(org_mcount[i] == m->count(i));

  // do intensive ghosting test
  for (int brg_dim=mesh_dim-1; brg_dim>=0; --brg_dim)
    for (int num_layer=1; num_layer<=3; ++num_layer)
      for (int include_copy=0; include_copy<=1; ++include_copy)
      {
        int before_mcount = m->count(mesh_dim);
        pumi_ghost_createlayer (m, brg_dim, mesh_dim, num_layer, include_copy);
        int total_mcount_diff=0, mcount_diff = m->count(mesh_dim)-before_mcount;
        MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm());
        if (!pumi_rank()) std::cout<<"\n[test_pumi] pumi_ghost_createlayer (bd "<<brg_dim<<", gd "<<mesh_dim<<", nl "<<num_layer<<", ic"<<include_copy<<") #ghost increase="<<total_mcount_diff<<"\n";
        pumi_mesh_print(m);
      }
  pumi_ghost_delete(m);

  for (int i=0; i<4; ++i)
  {
    if (org_mcount[i] != m->count(i)) 
       std::cout<<"("<<pumi_rank()<<") ERROR dim "<<i<<": org ent count "<<org_mcount[i]<<", current ent count "<<m->count(i)<<"\n";
//    assert(org_mcount[i] == m->count(i));
  }
  
  delete [] org_mcount;

  // print elapsed time and increased heap memory
  pumi_printtimemem("[test_pumi] elapsed time and increased heap memory:", pumi_gettime()-begin_time, pumi_getmem()-begin_mem);

  // clean-up
  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
}

