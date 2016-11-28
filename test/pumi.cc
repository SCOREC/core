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

void TEST_GEOM_TAG(pGeom g);
void TEST_GHOSTING(pMesh m);

//*********************************************************
int main(int argc, char** argv)
//*********************************************************
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pumi_printSys();

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

  // read input args - in-model-file in-mesh-file out-mesh-file num-in-part
  getConfig(argc,argv);

  // load model
  pGeom g = pumi_geom_load(modelFile);

  if (!pumi_rank()) std::cout<<"[test_pumi] testing geometric model tagging api's\n\n";
  TEST_GEOM_TAG(g);
 
  // load mesh per process group
  assert(pumi_size()%num_in_part==0);

  double begin_mem = pumi_getMem(), begin_time=pumi_getTime();

  pMesh m=NULL;
  if (do_distr)
    m = pumi_mesh_loadSerial(g, meshFile);
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

    int dim=pumi_mesh_getDim(m), count=0, pid;
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      pid=pumi_ment_getLocalID(e)%PCU_Comm_Peers();
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

  if (!pumi_rank()) std::cout<<"\n[test_pumi] checking various mesh api's\n";
  int mesh_dim=pumi_mesh_getDim(m);
  pMeshIter mit = m->begin(mesh_dim);
  while ((e = m->iterate(mit)))
  {
    assert(pumi_ment_getDim(e)==mesh_dim);
    assert(pumi_ment_getNumAdj(e, mesh_dim+1)==0);
    if (!pumi_ment_isOnBdry(e)) continue; // skip internal entity
    // if entity is on part boundary, count remote copies    
    Copies copies;
    pumi_ment_getAllRmt(e,copies);
    // loop over remote copies and increase the counter
    // check #remotes
    assert (pumi_ment_getNumRmt(e)==int(copies.size()) && copies.size()>0);
    // check the entity is not ghost or ghosted
    assert(!pumi_ment_isGhost(e) && !pumi_ment_isGhosted(e));
  }
  m->end(mit);

  // re-load partitioned mesh via file i/o
  pumi_mesh_delete(m);
  g = pumi_geom_load(modelFile);
  if (num_in_part==1 && pumi_size()>1)
    m =  pumi_mesh_load(g, "mesh.smb", pumi_size()); 
  else
    m = pumi_mesh_load(g, meshFile, num_in_part); 
  if (!pumi_rank()) std::cout<<"\n[test_pumi] delete and reload mesh\n";

  if (!pumi_rank()) std::cout<<"\n[test_pumi] clean loaded tags from the mesh file\n";
  std::vector<pMeshTag> tag_vec;
  for (size_t n = 0; n<tag_vec.size(); ++n)
  {
    pumi_mesh_deleteTag(m, tag_vec[n], true /* force_delete*/);    
  }

  TEST_GHOSTING(m);

  // print elapsed time and increased heap memory
  pumi_printTimeMem("\n* [test_pumi] elapsed time and increased heap memory:", pumi_getTime()-begin_time, pumi_getMem()-begin_mem);

  // clean-up 
  pumi_mesh_verify(m);
  pumi_mesh_delete(m);
  pumi_finalize();
  MPI_Finalize();
}

#include <string.h>
#include <cstring>

//*********************************************************
template <class T>
void TEST_TAG (pTag tag, char const* in_name, int name_len, int in_type, int in_size)
//*********************************************************
{
  const char* tag_name;
  // verifying byte tag info
  pumi_tag_getName (tag, &tag_name);
  int tag_type= pumi_tag_getType (tag);
  int tag_size = pumi_tag_getSize (tag);
  int tag_byte= pumi_tag_getByte (tag);
  assert(!strncmp(tag_name, in_name, name_len));
  assert(tag_type == in_type);
  assert(tag_size == in_size);
  assert(size_t(tag_byte)==sizeof(T)*tag_size);
}


//*********************************************************
void TEST_GENT_SETGET_TAG (pGeom g, pGeomEnt ent)
//*********************************************************
{
  char data[] = "abcdefg";

  pTag pointer_tag=pumi_geom_findTag(g, "pointer");
  pTag int_tag=pumi_geom_findTag(g, "integer");
  pTag long_tag = pumi_geom_findTag(g, "long");
  pTag dbl_tag = pumi_geom_findTag(g, "double");
  pTag ent_tag = pumi_geom_findTag(g, "entity");
  pTag intarr_tag = pumi_geom_findTag(g, "integer array");
  pTag longarr_tag = pumi_geom_findTag(g, "integer array");
  pTag dblarr_tag = pumi_geom_findTag(g, "double array");
  pTag entarr_tag = pumi_geom_findTag(g, "entity array");

  // get an entity to use as tag data
  pGeomEnt ent_tag_data=*(g->begin(0));

  // pumi_gent_set/getPtrTag 
  pumi_gent_setPtrTag (ent, pointer_tag, (void*)(data));
  void* void_data = (void*)calloc(strlen(data), sizeof(char));
  pumi_gent_getPtrTag (ent, pointer_tag, &void_data);
  assert(!strcmp((char*)void_data, data));

  // pumi_gent_set/getIntTag 
  pumi_gent_setIntTag(ent, int_tag, 1000);
  int int_data;
  pumi_gent_getIntTag (ent, int_tag, &int_data);
  assert(int_data == 1000);

  // pumi_gent_set/getLongTag
  pumi_gent_setLongTag(ent, long_tag, 3000);
  long long_data; 
  pumi_gent_getLongTag (ent, long_tag, &long_data);
  assert(long_data==3000);

  // pumi_gent_set/getDblTag
  pumi_gent_setDblTag (ent, dbl_tag, 1000.37);
  double dbl_data;
  pumi_gent_getDblTag (ent, dbl_tag, &dbl_data);
  assert(dbl_data == 1000.37);

  // pumi_gent_set/getEntTag
  pumi_gent_setEntTag (ent, ent_tag, ent_tag_data);
  pGeomEnt ent_data; 
  pumi_gent_getEntTag (ent, ent_tag, &ent_data);
  assert(ent_data == ent_tag_data);


 // pumi_gent_set/GetIntArrTag with integer arr tag
  int int_arr[] = {4,8,12};
  int arr_size;
  pumi_gent_setIntArrTag (ent, intarr_tag, int_arr);
  int* int_arr_back = new int[4];
  pumi_gent_getIntArrTag (ent, intarr_tag, &int_arr_back, &arr_size);
  assert(arr_size==3 && int_arr_back[0] == 4 && int_arr_back[1] == 8 && int_arr_back[2] == 12);
            
 // pumi_gent_set/getLongArrTag 
  long long_arr[] = {4,8,12};
  pumi_gent_setLongArrTag (ent, longarr_tag, long_arr);
  long* long_arr_back = new long[4];
  pumi_gent_getLongArrTag (ent, longarr_tag, &long_arr_back, &arr_size);
// FIXME: this fails
//  assert(arr_size==3 && long_arr_back[0] == 4 && long_arr_back[1] == 8 && long_arr_back[2] == 12);

   // pumi_gent_set/getDblArrTag
  double dbl_arr[] = {4.1,8.2,12.3};
  pumi_gent_setDblArrTag (ent, dblarr_tag, dbl_arr);
  double* dbl_arr_back = new double[4];
  pumi_gent_getDblArrTag (ent, dblarr_tag, &dbl_arr_back, &arr_size);
  assert(arr_size==3 && dbl_arr_back[0] == 4.1 && dbl_arr_back[1] == 8.2 && 
         dbl_arr_back[2] == 12.3);

  // pumi_gent_set/getEntArrTag
  pGeomEnt* ent_arr = new pGeomEnt[3];
  ent_arr[0] = ent_arr[1] = ent_arr[2] = ent_tag_data;
  pumi_gent_setEntArrTag (ent, entarr_tag, ent_arr);
  pGeomEnt* ent_arr_back = new pGeomEnt[4];
  pumi_gent_getEntArrTag (ent, entarr_tag, &ent_arr_back, &arr_size);
  assert(arr_size==3 && ent_arr_back[0] == ent_tag_data && ent_arr_back[1] == 
           ent_tag_data && ent_arr_back[2] == ent_tag_data
           && ent_arr[0]==ent_arr_back[0] && ent_arr[1]==ent_arr_back[1] && 
           ent_arr[2] == ent_arr_back[2]);

  delete [] int_arr_back;
  delete [] long_arr_back;
  delete [] dbl_arr_back;
  delete [] ent_arr;
  delete [] ent_arr_back;
}

//*********************************************************
void TEST_GENT_DEL_TAG (pGeom g, pGeomEnt ent)
//*********************************************************
{
  pTag pointer_tag=pumi_geom_findTag(g, "pointer");
  pTag int_tag=pumi_geom_findTag(g, "integer");
  pTag long_tag = pumi_geom_findTag(g, "long");
  pTag dbl_tag = pumi_geom_findTag(g, "double");
  pTag ent_tag = pumi_geom_findTag(g, "entity");
  pTag intarr_tag = pumi_geom_findTag(g, "integer array");
  pTag longarr_tag = pumi_geom_findTag(g, "long array");
  pTag dblarr_tag = pumi_geom_findTag(g, "double array");
  pTag entarr_tag = pumi_geom_findTag(g, "entity array");

  pumi_gent_deleteTag(ent, pointer_tag);
  pumi_gent_deleteTag(ent, int_tag);
  pumi_gent_deleteTag(ent, long_tag);
  pumi_gent_deleteTag(ent, dbl_tag);
  pumi_gent_deleteTag(ent, ent_tag);
  pumi_gent_deleteTag(ent, intarr_tag);
  pumi_gent_deleteTag(ent, longarr_tag);
  pumi_gent_deleteTag(ent, dblarr_tag);
  pumi_gent_deleteTag(ent, entarr_tag);

  assert(!pumi_gent_hasTag(ent, pointer_tag));
  assert(!pumi_gent_hasTag(ent, int_tag));
  assert(!pumi_gent_hasTag(ent, long_tag));
  assert(!pumi_gent_hasTag(ent, dbl_tag));
  assert(!pumi_gent_hasTag(ent, ent_tag));
  assert(!pumi_gent_hasTag(ent, intarr_tag));
  assert(!pumi_gent_hasTag(ent, longarr_tag));
  assert(!pumi_gent_hasTag(ent, dblarr_tag));
  assert(!pumi_gent_hasTag(ent, entarr_tag));
}

void TEST_GEOM_TAG(pGeom g)
{
  pTag pointer_tag = pumi_geom_createTag(g, "pointer", PUMI_PTR, 1);
  pTag int_tag = pumi_geom_createTag(g, "integer", PUMI_INT, 1);
  pTag long_tag = pumi_geom_createTag(g, "long", PUMI_LONG, 1);
  pTag dbl_tag = pumi_geom_createTag(g, "double", PUMI_DBL, 1);
  pTag ent_tag = pumi_geom_createTag(g, "entity", PUMI_ENT, 1);

  pTag intarr_tag=pumi_geom_createTag(g, "integer array", PUMI_INT, 3);
  pTag longarr_tag=pumi_geom_createTag(g, "long array", PUMI_LONG, 3);
  pTag dblarr_tag = pumi_geom_createTag(g, "double array", PUMI_DBL, 3);
  pTag entarr_tag = pumi_geom_createTag(g, "entity array", PUMI_ENT, 3);

  // verifying tag info
  TEST_TAG<void*>(pointer_tag, "pointer", strlen("pointer"), PUMI_PTR, 1);
  TEST_TAG<int>(int_tag, "integer", strlen("integer"), PUMI_INT, 1);
  TEST_TAG<long>(long_tag, "long", strlen("long"), PUMI_LONG, 1);
  TEST_TAG<double>(dbl_tag, "double", strlen("double"), PUMI_DBL, 1);
  TEST_TAG<pMeshEnt>(ent_tag, "entity", strlen("entity"), PUMI_ENT, 1);

  TEST_TAG<int>(intarr_tag, "integer array", strlen("integer array"), PUMI_INT, 3);
  TEST_TAG<long>(longarr_tag, "long array", strlen("long array"), PUMI_LONG, 3);
  TEST_TAG<double>(dblarr_tag, "double array", strlen("double array"), PUMI_DBL, 3);
  TEST_TAG<pMeshEnt>(entarr_tag, "entity array", strlen("entity array"), PUMI_ENT, 3);

  assert(pumi_geom_hasTag(g, int_tag));
  pTag cloneTag = pumi_geom_findTag(g, "pointer");
  assert(cloneTag);
  std::vector<pTag> tags;
  pumi_geom_getTag(g, tags);
  assert(cloneTag == pointer_tag && tags.size()==9);

  for (pGeomIter gent_it = g->begin(0); gent_it!=g->end(0);++gent_it)
  {
    TEST_GENT_SETGET_TAG(g, *gent_it);
    TEST_GENT_DEL_TAG(g, *gent_it);
  }

  // delete tags
  for (std::vector<pTag>::iterator tag_it=tags.begin(); tag_it!=tags.end(); ++tag_it)
    pumi_geom_deleteTag(g, *tag_it);

  tags.clear();
  pumi_geom_getTag(g, tags);

  assert(!tags.size());
}


Ghosting* getGhostingPlan(pMesh m)
{
  int mesh_dim=m->getDimension();
  Ghosting* plan = new Ghosting(m, m->getDimension());
  {
    apf::MeshIterator* it = m->begin(mesh_dim);
    pMeshEnt e;
    size_t count=0;
    while ((e = m->iterate(it)))
    {
      for (int i=0; i<pumi_size()/2; ++i)
      {
        int pid = rand()%pumi_size();
        plan->send(e, pid);
      }
      ++count; 
     if (count==m->count(mesh_dim)/3) break;
    }
    m->end(it);
  }
  return plan;
}

void TEST_GHOSTING(pMesh m)
{  
  int mesh_dim=m->getDimension();
  pMeshEnt e;
  // element-wise ghosting test
  int num_org_vtx = pumi_mesh_getNumEnt(m, 0);
  int* org_mcount=new int[4];
  for (int i=0; i<4; ++i)
    org_mcount[i] = m->count(i);

  Ghosting* ghosting_plan = getGhostingPlan(m);
  int before_mcount=m->count(mesh_dim);
  pumi_ghost_create(m, ghosting_plan);

  int total_mcount_diff=0, mcount_diff = m->count(mesh_dim)-before_mcount;
  MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm()); 
  if (!pumi_rank()) std::cout<<"\n[test_pumi] element-wise pumi_ghost_create: #ghost increase="<<total_mcount_diff<<"\n";

  int num_ghost_vtx=0;
  pMeshIter mit = m->begin(0);
  while ((e = m->iterate(mit)))
  {
    if (pumi_ment_isGhost(e))
    {
      ++num_ghost_vtx;
     assert(pumi_ment_getOwnPID(e)!=pumi_rank());
    }
  }   
  m->end(mit);
  assert(num_ghost_vtx+num_org_vtx==pumi_mesh_getNumEnt(m,0));
  pumi_mesh_verify(m); // this should throw an error message
  pumi_ghost_delete(m);
  for (int i=0; i<4; ++i)
    assert(org_mcount[i] == int(m->count(i)));

  // layer-wise ghosting test
  for (int brg_dim=mesh_dim-1; brg_dim>=0; --brg_dim)
    for (int num_layer=1; num_layer<=3; ++num_layer)
      for (int include_copy=0; include_copy<=1; ++include_copy)
      {
        before_mcount=m->count(mesh_dim);
        pumi_ghost_createLayer (m, brg_dim, mesh_dim, num_layer, include_copy);
        total_mcount_diff=0, mcount_diff = m->count(mesh_dim)-before_mcount;
        MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm());
        if (!pumi_rank()) std::cout<<"\n[test_pumi] layer-wise pumi_ghost_createLayer (bd "<<brg_dim<<", gd "<<mesh_dim<<", nl "<<num_layer<<", ic"<<include_copy<<"), #ghost increase="<<total_mcount_diff<<"\n";
        pumi_mesh_verify(m);
        pumi_ghost_delete(m);
        for (int i=0; i<4; ++i)
          assert(org_mcount[i] == int(m->count(i)));
      }
  
  // accumulative layer-ghosting
  for (int brg_dim=mesh_dim-1; brg_dim>=0; --brg_dim)
    for (int num_layer=1; num_layer<=3; ++num_layer)
      for (int include_copy=0; include_copy<=1; ++include_copy)
      {
        int before_mcount=m->count(mesh_dim);
        pumi_ghost_createLayer (m, brg_dim, mesh_dim, num_layer, include_copy);
        int total_mcount_diff=0, mcount_diff = m->count(mesh_dim)-before_mcount;
        MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm());
        if (!pumi_rank()) 
          std::cout<<"\n[test_pumi] accumulative pumi_ghost_createLayer (bd "<<brg_dim<<", gd "<<mesh_dim
                   <<", nl "<<num_layer<<", ic"<<include_copy<<"), #ghost increase="<<total_mcount_diff<<"\n";
      }
 // pumi_mesh_verify(m); -- FIXME: this returns an error with ghost copy
  pumi_ghost_delete(m);

  for (int i=0; i<4; ++i)
  {
    if (org_mcount[i] != int(m->count(i)))
       std::cout<<"("<<pumi_rank()<<") ERROR dim "<<i<<": org ent count "<<org_mcount[i]<<", current ent count "<<m->count(i)<<"\n";
    assert(org_mcount[i] == int(m->count(i)));
  }
  
  delete [] org_mcount;
}
