#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfZoltan.h>
#include <pcu_util.h>
#include <cstdlib>
#include <pumi.h>
#include <unistd.h>
#include <iostream>
#include <algorithm>
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

bool is_double_isequal(double A, double B)
{
  double maxDiff = 1e-5;
  double maxRelDiff = 1e-5;
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ 
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;
 
    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}


void TEST_GEOM_TAG(pGeom g);
void TEST_MESH(pMesh m);
void TEST_MESH_TAG(pMesh m);
void TEST_NEW_MESH(pMesh m);
void TEST_GHOSTING(pMesh m);
void TEST_FIELD(pMesh m);

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

  if (!pumi_rank()) std::cout<<"[test_pumi] testing geometric model/entity api's\n\n";
  {
    // test geom_find and gent_adj
    for (pGeomIter gent_it = g->begin(2); gent_it!=g->end(2);++gent_it)
    {
      int id = pumi_gent_getID(*gent_it);
      pGeomEnt ge = pumi_geom_findEnt(g, 2, id);
      PCU_ALWAYS_ASSERT(ge == *gent_it);
      std::vector<pGeomEnt> adj_edges;
      pumi_gent_getAdj(ge, 1, adj_edges);
      std::vector<pGeomEnt> adj_vertices;
      pumi_gent_getAdj(ge, 0, adj_vertices);

      for (size_t nv=0; nv<adj_vertices.size(); ++nv)
      {
        pGeomEnt gv=adj_vertices.at(nv);
        std::vector<pGeomEnt> adj_faces;
        pumi_gent_getAdj(gv, 2, adj_faces);
        PCU_ALWAYS_ASSERT(std::find(adj_faces.begin(), adj_faces.end(), ge)!=adj_faces.end());
      }

    }

    TEST_GEOM_TAG(g);
  }
 
  // load mesh per process group
  PCU_ALWAYS_ASSERT(pumi_size()%num_in_part==0);

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
    pMeshIter it = m->begin(dim);
    while ((e = m->iterate(it)))
    {
      pid=pumi_ment_getID(e)%PCU_Comm_Peers();
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

  TEST_MESH(m);
 
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
    pumi_mesh_deleteTag(m, tag_vec[n], true /* force_delete*/);    

  int num_mesh_rg=pumi_mesh_getNumEnt(m,3);
  PCU_ALWAYS_ASSERT(pumi_mesh_findEnt(m, 3, num_mesh_rg+1)==NULL);

  // create global ID
  pumi_mesh_createGlobalID(m);

  TEST_FIELD(m);

  int num_field=pumi_mesh_getNumField(m);
  std::vector<pField> fields;
  for (int i=0; i<num_field;++i)
    fields.push_back(pumi_mesh_getField(m, i));

  for (std::vector<pField>::iterator fit=fields.begin(); fit!=fields.end(); ++fit)
    pumi_field_freeze(*fit);

  if (!pumi_rank()) std::cout<<"\n[test_pumi] "<<fields.size()<<" field(s) generated, synchronized, and frozen\n\n";

  TEST_GHOSTING(m);

  // delete numbering and ID
  std::vector<pGlobalNumbering> numberings;
  int num_gn=pumi_mesh_getNumGlobalNumbering(m);
  for (int i=0; i<num_gn;++i)
    numberings.push_back(pumi_mesh_getGlobalNumbering(m, i));
  PCU_ALWAYS_ASSERT(pumi_mesh_getNumGlobalNumbering(m)==(int)numberings.size());

  for (int i=0; i<(int)numberings.size(); ++i)
    pumi_numbering_deleteGlobal(numberings.at(i));
  pumi_mesh_deleteGlobalID(m);

  // delete fields
  // FIXME: FieldShape doesn't get removed along the field
  for (std::vector<pField>::iterator fit=fields.begin(); fit!=fields.end(); ++fit)
    pumi_field_delete(*fit);
  if (!pumi_rank()) std::cout<<"\n[test_pumi] field and numbering deleted\n";
  pumi_mesh_verify(m);

  TEST_MESH_TAG(m);

  // clean-up 
  pumi_mesh_delete(m);

  // print elapsed time and increased heap memory
  pumi_printTimeMem("\n* [test_pumi] elapsed time and increased heap memory:", pumi_getTime()-begin_time, pumi_getMem()-begin_mem);

  pumi_finalize();
  MPI_Finalize();
}

void TEST_MESH(pMesh m)
{
  int mesh_dim=pumi_mesh_getDim(m);
  pMeshEnt e;
  std::vector<pMeshEnt> adj_vtx;
  std::vector<pMeshEnt> adj_elem;

  pMeshIter mit = m->begin(mesh_dim);
  while ((e = m->iterate(mit)))
  {
    PCU_ALWAYS_ASSERT(pumi_ment_getDim(e)==mesh_dim);
    PCU_ALWAYS_ASSERT(pumi_ment_getNumAdj(e, mesh_dim+1)==0);
    // check adjacency
    adj_vtx.clear();
    pumi_ment_getAdj(e, 0, adj_vtx);
    PCU_ALWAYS_ASSERT((size_t)pumi_ment_getNumAdj(e, 0)==adj_vtx.size());

    adj_elem.clear();
    pumi_ment_getAdj(adj_vtx.at(0), mesh_dim, adj_elem);
    PCU_ALWAYS_ASSERT(std::find(adj_elem.begin(), adj_elem.end(), e)!=adj_elem.end());

    apf::Adjacent adjacent;
    m->getAdjacent(adj_vtx.at(0), mesh_dim,adjacent);      
    PCU_ALWAYS_ASSERT(adjacent.getSize()==adj_elem.size());

    if (!pumi_ment_isOnBdry(e)) continue; // skip internal entity
    // if entity is on part boundary, count remote copies    
    Copies copies;
    pumi_ment_getAllRmt(e,copies);
    // loop over remote copies and increase the counter
    // check #remotes
    PCU_ALWAYS_ASSERT(pumi_ment_getNumRmt(e)==int(copies.size()) && copies.size()>0);
    Parts parts;
    pumi_ment_getResidence(e, parts);
    PCU_ALWAYS_ASSERT(parts.size()==copies.size()+1);

    for (pCopyIter it = copies.begin();
           it != copies.end(); ++it)
      PCU_ALWAYS_ASSERT(it->first!=pumi_rank());
    // check the entity is not ghost or ghosted
    PCU_ALWAYS_ASSERT(!pumi_ment_isGhost(e) && !pumi_ment_isGhosted(e));
  }
  m->end(mit);

  if (!pumi_rank()) std::cout<<"\n[test_pumi] testing  mesh/entity apis\n";
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
  PCU_ALWAYS_ASSERT(!strncmp(tag_name, in_name, name_len));
  PCU_ALWAYS_ASSERT(tag_type == in_type);
  PCU_ALWAYS_ASSERT(tag_size == in_size);
  PCU_ALWAYS_ASSERT(size_t(tag_byte)==sizeof(T)*tag_size);
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
  PCU_ALWAYS_ASSERT(!strcmp((char*)void_data, data));

  // pumi_gent_set/getIntTag 
  pumi_gent_setIntTag(ent, int_tag, 1000);
  int int_data;
  pumi_gent_getIntTag (ent, int_tag, &int_data);
  PCU_ALWAYS_ASSERT(int_data == 1000);

  // pumi_gent_set/getLongTag
  pumi_gent_setLongTag(ent, long_tag, 3000);
  long long_data; 
  pumi_gent_getLongTag (ent, long_tag, &long_data);
  PCU_ALWAYS_ASSERT(long_data==3000);

  // pumi_gent_set/getDblTag
  pumi_gent_setDblTag (ent, dbl_tag, 1000.37);
  double dbl_data;
  pumi_gent_getDblTag (ent, dbl_tag, &dbl_data);
  PCU_ALWAYS_ASSERT(is_double_isequal(dbl_data,1000.37));

  // pumi_gent_set/getEntTag
  pumi_gent_setEntTag (ent, ent_tag, ent_tag_data);
  pGeomEnt ent_data; 
  pumi_gent_getEntTag (ent, ent_tag, &ent_data);
  PCU_ALWAYS_ASSERT(ent_data == ent_tag_data);


 // pumi_gent_set/GetIntArrTag with integer arr tag
  int int_arr[] = {4,8,12};
  int arr_size;
  pumi_gent_setIntArrTag (ent, intarr_tag, int_arr);
  int* int_arr_back = new int[4];
  pumi_gent_getIntArrTag (ent, intarr_tag, &int_arr_back, &arr_size);
  PCU_ALWAYS_ASSERT(arr_size==3 && int_arr_back[0] == int_arr[0] && int_arr_back[1] == int_arr[1] && int_arr_back[2] == int_arr[2]);
            
 // pumi_gent_set/getLongArrTag 
  long long_arr[] = {4,8,12};
  pumi_gent_setLongArrTag (ent, longarr_tag, long_arr);
  long* long_arr_back = new long[4];
  pumi_gent_getLongArrTag (ent, longarr_tag, &long_arr_back, &arr_size);
// FIXME: this fails
//  PCU_ALWAYS_ASSERT(arr_size==3 && long_arr_back[0] == long_arr[0] 
//         && long_arr_back[1] == long_arr[1] && long_arr_back[2] == long_arr[2]);

   // pumi_gent_set/getDblArrTag
  double dbl_arr[] = {4.1,8.2,12.3};
  pumi_gent_setDblArrTag (ent, dblarr_tag, dbl_arr);
  double* dbl_arr_back = new double[4];
  pumi_gent_getDblArrTag (ent, dblarr_tag, &dbl_arr_back, &arr_size);
  PCU_ALWAYS_ASSERT(arr_size==3 && dbl_arr_back[0] == dbl_arr[0] && dbl_arr_back[1] == dbl_arr[1] && 
         dbl_arr_back[2] == dbl_arr[2]);

  // pumi_gent_set/getEntArrTag
  pGeomEnt* ent_arr = new pGeomEnt[3];
  ent_arr[0] = ent_arr[1] = ent_arr[2] = ent_tag_data;
  pumi_gent_setEntArrTag (ent, entarr_tag, ent_arr);
  pGeomEnt* ent_arr_back = new pGeomEnt[4];
  pumi_gent_getEntArrTag (ent, entarr_tag, &ent_arr_back, &arr_size);
  PCU_ALWAYS_ASSERT(arr_size==3 && ent_arr_back[0] == ent_tag_data && ent_arr_back[1] == 
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

  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, pointer_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, int_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, long_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, dbl_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, ent_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, intarr_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, longarr_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, dblarr_tag));
  PCU_ALWAYS_ASSERT(!pumi_gent_hasTag(ent, entarr_tag));
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

  PCU_ALWAYS_ASSERT(pumi_geom_hasTag(g, int_tag));
  pTag cloneTag = pumi_geom_findTag(g, "pointer");
  PCU_ALWAYS_ASSERT(cloneTag);
  std::vector<pTag> tags;
  pumi_geom_getTag(g, tags);
  PCU_ALWAYS_ASSERT(cloneTag == pointer_tag && tags.size()==9);

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

  PCU_ALWAYS_ASSERT(!tags.size());
}

//*********************************************************
void TEST_MENT_SETGET_TAG (pMesh m, pMeshEnt ent)
//*********************************************************
{
  pMeshTag int_tag=pumi_mesh_findTag(m, "integer");
  pMeshTag long_tag = pumi_mesh_findTag(m, "long");
  pMeshTag dbl_tag = pumi_mesh_findTag(m, "double");
  pMeshTag intarr_tag = pumi_mesh_findTag(m, "integer array");
  pMeshTag longarr_tag = pumi_mesh_findTag(m, "long array");
  pMeshTag dblarr_tag = pumi_mesh_findTag(m, "double array");

  // pumi_ment_set/getIntTag 
  int int_value=pumi_ment_getID(ent), int_data;
  pumi_ment_setIntTag(ent, int_tag, &int_value);
  pumi_ment_getIntTag (ent, int_tag, &int_data);
  PCU_ALWAYS_ASSERT(int_data == int_value);

  // pumi_ment_set/getLongTag
  long long_value=3000, long_data;
  pumi_ment_setLongTag(ent, long_tag, &long_value);
  pumi_ment_getLongTag (ent, long_tag, &long_data);
  PCU_ALWAYS_ASSERT(long_data==long_value);

  // pumi_gent_set/getDblTag
  double dbl_value=1000.37, dbl_data;
  pumi_ment_setDblTag (ent, dbl_tag, &dbl_value);
  pumi_ment_getDblTag (ent, dbl_tag, &dbl_data);
  PCU_ALWAYS_ASSERT(is_double_isequal(dbl_data,dbl_value));

 // pumi_ment_set/getIntTag with integer arr
  int int_arr[] = {4,8,12};
  pumi_ment_setIntTag (ent, intarr_tag, int_arr);
  int* int_arr_back = new int[3];
  pumi_ment_getIntTag (ent, intarr_tag, int_arr_back);
  PCU_ALWAYS_ASSERT(int_arr_back[0] == 4 && int_arr_back[1] == 8 && int_arr_back[2] == 12);
            
 // pumi_ment_set/getLongTag with long arr
  long long_arr[] = {4,8,12};
  pumi_ment_setLongTag (ent, longarr_tag, long_arr);
  long* long_arr_back = new long[3];
  pumi_ment_getLongTag (ent, longarr_tag, long_arr_back);
  PCU_ALWAYS_ASSERT(long_arr_back[0] == long_arr[0] && long_arr_back[1] == long_arr[1] && long_arr_back[2] == long_arr[2]);

   // pumi_ment_set/getDblTag with double arr
  double dbl_arr[] = {4.1,8.2,12.3};
  pumi_ment_setDblTag (ent, dblarr_tag, dbl_arr);
  double* dbl_arr_back = new double[3];
  pumi_ment_getDblTag (ent, dblarr_tag, dbl_arr_back);
  PCU_ALWAYS_ASSERT(dbl_arr_back[0] == dbl_arr[0] && dbl_arr_back[1] == dbl_arr[1] && dbl_arr_back[2] == dbl_arr[2]);

  delete [] int_arr_back;
  delete [] long_arr_back;
  delete [] dbl_arr_back;
}

//*********************************************************
void TEST_MENT_DEL_TAG (pMesh m, pMeshEnt ent)
//*********************************************************
{
  pMeshTag int_tag=pumi_mesh_findTag(m, "integer");
  pMeshTag long_tag = pumi_mesh_findTag(m, "long");
  pMeshTag dbl_tag = pumi_mesh_findTag(m, "double");
  pMeshTag intarr_tag = pumi_mesh_findTag(m, "integer array");
  pMeshTag longarr_tag = pumi_mesh_findTag(m, "long array");
  pMeshTag dblarr_tag = pumi_mesh_findTag(m, "double array");

  pumi_ment_deleteTag(ent, int_tag);
  pumi_ment_deleteTag(ent, long_tag);
  pumi_ment_deleteTag(ent, dbl_tag);
  pumi_ment_deleteTag(ent, intarr_tag);
  pumi_ment_deleteTag(ent, longarr_tag);
  pumi_ment_deleteTag(ent, dblarr_tag);

  PCU_ALWAYS_ASSERT(!pumi_ment_hasTag(ent, int_tag));
  PCU_ALWAYS_ASSERT(!pumi_ment_hasTag(ent, long_tag));
  PCU_ALWAYS_ASSERT(!pumi_ment_hasTag(ent, dbl_tag));
  PCU_ALWAYS_ASSERT(!pumi_ment_hasTag(ent, intarr_tag));
  PCU_ALWAYS_ASSERT(!pumi_ment_hasTag(ent, longarr_tag));
  PCU_ALWAYS_ASSERT(!pumi_ment_hasTag(ent, dblarr_tag));
}

//*********************************************************
void TEST_MESH_TAG(pMesh m)
//*********************************************************
{
  pMeshTag int_tag = pumi_mesh_createIntTag(m, "integer", 1);
  pMeshTag long_tag = pumi_mesh_createLongTag(m, "long", 1);
  pMeshTag dbl_tag = pumi_mesh_createDblTag(m, "double", 1);

  pMeshTag intarr_tag=pumi_mesh_createIntTag(m, "integer array", 3);
  pMeshTag longarr_tag=pumi_mesh_createLongTag(m, "long array", 3);
  pMeshTag dblarr_tag = pumi_mesh_createDblTag(m, "double array", 3);

  PCU_ALWAYS_ASSERT(pumi_mesh_hasTag(m, int_tag));
  PCU_ALWAYS_ASSERT(pumi_mesh_hasTag(m, long_tag));
  PCU_ALWAYS_ASSERT(pumi_mesh_hasTag(m, dbl_tag));

  PCU_ALWAYS_ASSERT(pumi_mesh_hasTag(m, intarr_tag));
  PCU_ALWAYS_ASSERT(pumi_mesh_hasTag(m, longarr_tag));
  PCU_ALWAYS_ASSERT(pumi_mesh_hasTag(m, dblarr_tag));

  pMeshTag cloneTag = pumi_mesh_findTag(m, "double");
  PCU_ALWAYS_ASSERT(cloneTag);
  std::vector<pMeshTag> tags;
  pumi_mesh_getTag(m, tags);
  PCU_ALWAYS_ASSERT(cloneTag == dbl_tag);

  pMeshIter it = m->begin(0);
  pMeshEnt e;

  while ((e = m->iterate(it)))
  {
    TEST_MENT_SETGET_TAG(m, e);
    TEST_MENT_DEL_TAG(m, e);
  }
  m->end(it);

  // delete tags
  for (std::vector<pMeshTag>::iterator tag_it=tags.begin(); tag_it!=tags.end(); ++tag_it)
    pumi_mesh_deleteTag(m, *tag_it);

  tags.clear();
  pumi_mesh_getTag(m, tags);

  PCU_ALWAYS_ASSERT(!tags.size());
}


void TEST_NEW_MESH(pMesh m)
{
  // change chape of the mesh
  pShape s = pumi_mesh_getShape(m);
  PCU_ALWAYS_ASSERT(pumi_shape_getNumNode(s, 1)==0);

  pumi_mesh_setShape(m, pumi_shape_getLagrange(2));
  PCU_ALWAYS_ASSERT(pumi_shape_getNumNode(pumi_mesh_getShape(m), 1)==1);

  // create an empty mesh
  pGeom new_g = pumi_geom_load (NULL, "null");
  pMesh new_m = pumi_mesh_create(new_g, 2);

  double xyz[3];
  pMeshIter it = m->begin(1);
  pMeshEnt e;
  std::vector<pMeshEnt> new_vertices;
  std::vector<pMeshEnt> new_edges;

  while ((e = m->iterate(it)))
  {
    pumi_node_getCoord(e, 0, xyz);
    new_vertices.push_back (pumi_mesh_createVtx(new_m, NULL, xyz));
  }
  m->end(it);

  for (size_t i=0; i< new_vertices.size()/2-1; ++i)
  {
    pMeshEnt vertices[2];
    vertices[0] = new_vertices[i];    
    vertices[1] = new_vertices[i+1]; 
    new_edges.push_back(pumi_mesh_createEnt(new_m, NULL, 1, vertices));
  }

  // note we ignore face-edge order in this example
  for (size_t i=0; i< new_edges.size()/3-1; ++i)
  {
    pMeshEnt edges[3];
    edges[0] = new_edges[i];     
    edges[1] = new_edges[i+1]; 
    edges[2] = new_edges[i+2];
    pumi_mesh_createEnt(new_m, NULL, 2, edges);
  }

  for (int d=0; d<=pumi_mesh_getDim(new_m);++d)
    PCU_ALWAYS_ASSERT(pumi_mesh_getNumEnt(new_m,d));

//  pumi_mesh_freeze(new_m);

  if (!pumi_rank()) std::cout<<"\n[test_pumi] new mesh constructed (#v "
         <<pumi_mesh_getNumEnt(new_m, 0)
         <<", #e "<<pumi_mesh_getNumEnt(new_m, 1)
         <<", #f "<<pumi_mesh_getNumEnt(new_m, 2)
         <<", #r "<<pumi_mesh_getNumEnt(new_m, 3)<<")\n";

//  pumi_mesh_delete(new_m);
}

void TEST_FIELD(pMesh m)
{
  int num_dofs_per_node=3;
  pField f = pumi_mesh_findField(m, "xyz_field");
  pMeshIter it;
  pMeshEnt e;
  double data[3];
  double xyz[3];

  // create field and set the field data
  if (!f)
  {
    f=pumi_field_create(m, "xyz_field", num_dofs_per_node);
    // create global numbering
    pumi_numbering_createGlobal(m, "xyz_numbering", pumi_field_getShape(f));

    PCU_ALWAYS_ASSERT(pumi_field_getName(f)==std::string("xyz_field"));
    PCU_ALWAYS_ASSERT(pumi_field_getType(f)==PUMI_PACKED);
    PCU_ALWAYS_ASSERT(pumi_field_getSize(f)==num_dofs_per_node);
    it = m->begin(0);
    
    while ((e = m->iterate(it)))
    {
      if (!pumi_ment_isOwned(e)) continue;
      if (pumi_ment_isOnBdry(e)) 
        for (int i=0; i<3;++i) 
          xyz[i] = pumi_ment_getGlobalID(e);
      else 
        pumi_node_getCoord(e, 0, xyz);
      pumi_ment_setField(e, f, 0, xyz);
    }
    m->end(it);

    pumi_field_synchronize(f);
  }

  it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    pumi_node_getCoord(e, 0, xyz);
    pumi_ment_getField(e, f, 0, data);
    for (int i=0; i<3;++i) 
      if (pumi_ment_isOnBdry(e)) 
        PCU_ALWAYS_ASSERT(data[i] == pumi_ment_getGlobalID(e));
      else
        PCU_ALWAYS_ASSERT(data[i] == xyz[i]);
  }
  m->end(it);
  pumi_field_verify(m, f);
}

Ghosting* getGhostingPlan(pMesh m)
{
  int mesh_dim=pumi_mesh_getDim(m);
  Ghosting* plan = new Ghosting(m, pumi_mesh_getDim(m));
  {
    pMeshIter it = m->begin(mesh_dim);
    pMeshEnt e;
    int count=0;
    while ((e = m->iterate(it)))
    {
      for (int i=0; i<pumi_size()/2; ++i)
      {
        int pid = rand()%pumi_size();
        plan->send(e, pid);
      }
      ++count; 
     if (count==pumi_mesh_getNumEnt(m, mesh_dim)/3) break;
    }
    m->end(it);
  }
  return plan;
}

void TEST_GHOSTING(pMesh m)
{  
  int mesh_dim=pumi_mesh_getDim(m);
  pMeshEnt e;
  // element-wise ghosting test
  int num_org_vtx = pumi_mesh_getNumEnt(m, 0);
  int* org_mcount=new int[4];
  for (int i=0; i<4; ++i)
    org_mcount[i] = pumi_mesh_getNumEnt(m, i);
  
  Ghosting* ghosting_plan = getGhostingPlan(m);
  int before_mcount=pumi_mesh_getNumEnt(m, mesh_dim);
  pumi_ghost_create(m, ghosting_plan);

  int total_mcount_diff=0, mcount_diff = pumi_mesh_getNumEnt(m, mesh_dim)-before_mcount;
  MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm()); 
  if (!pumi_rank()) std::cout<<"\n[test_pumi] element-wise pumi_ghost_create: #ghost increase="<<total_mcount_diff<<"\n";

  int num_ghost_vtx=0;
  pMeshIter mit = m->begin(0);
  while ((e = m->iterate(mit)))
  {
    // FIXME: test pumi_ment_isOn and pumi_ment_getResidence with ghost copies
    if (pumi_ment_isGhost(e))
    {
      ++num_ghost_vtx;
     PCU_ALWAYS_ASSERT(pumi_ment_getOwnPID(e)!=pumi_rank());
    }
  }   
  m->end(mit);
  PCU_ALWAYS_ASSERT(num_ghost_vtx+num_org_vtx==pumi_mesh_getNumEnt(m,0));
  pumi_mesh_verify(m);
  TEST_FIELD(m);
  pumi_ghost_delete(m);
  
  for (int i=0; i<4; ++i)
    PCU_ALWAYS_ASSERT(org_mcount[i] == pumi_mesh_getNumEnt(m, i));

  // layer-wise ghosting test
  for (int brg_dim=mesh_dim-1; brg_dim>=0; --brg_dim)
    for (int num_layer=1; num_layer<=3; ++num_layer)
      for (int include_copy=0; include_copy<=1; ++include_copy)
      {
        before_mcount=pumi_mesh_getNumEnt(m, mesh_dim);
        pumi_ghost_createLayer (m, brg_dim, mesh_dim, num_layer, include_copy);
        total_mcount_diff=0, mcount_diff = pumi_mesh_getNumEnt(m, mesh_dim)-before_mcount;
        MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm());
        if (!pumi_rank()) std::cout<<"\n[test_pumi] layer-wise pumi_ghost_createLayer (bd "<<brg_dim<<", gd "<<mesh_dim<<", nl "<<num_layer<<", ic"<<include_copy<<"), #ghost increase="<<total_mcount_diff<<"\n";
        pumi_mesh_verify(m);
        TEST_FIELD(m);
        pumi_ghost_delete(m);
        for (int i=0; i<4; ++i)
          PCU_ALWAYS_ASSERT(org_mcount[i] == pumi_mesh_getNumEnt(m, i));
      }
  
  // accumulative layer-ghosting
  for (int brg_dim=mesh_dim-1; brg_dim>=0; --brg_dim)
    for (int num_layer=1; num_layer<=3; ++num_layer)
      for (int include_copy=0; include_copy<=1; ++include_copy)
      {
        int before_mcount=pumi_mesh_getNumEnt(m, mesh_dim);
        pumi_ghost_createLayer (m, brg_dim, mesh_dim, num_layer, include_copy);
        int total_mcount_diff=0, mcount_diff = pumi_mesh_getNumEnt(m, mesh_dim)-before_mcount;
        MPI_Allreduce(&mcount_diff, &total_mcount_diff,1, MPI_INT, MPI_SUM, PCU_Get_Comm());
        if (!pumi_rank()) 
          std::cout<<"\n[test_pumi] accumulative pumi_ghost_createLayer (bd "<<brg_dim<<", gd "<<mesh_dim
                   <<", nl "<<num_layer<<", ic"<<include_copy<<"), #ghost increase="<<total_mcount_diff<<"\n";
      }

  // test field accumulation with ghosted mesh  
  pField f = pumi_mesh_findField(m, "xyz_field");
  double data[3];
  double xyz[3];

  pumi_field_accumulate(f);

  pMeshIter it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    pumi_node_getCoord(e, 0, xyz);
    pumi_ment_getField(e, f, 0, data);
    for (int i=0; i<3;++i) 
      if (pumi_ment_isOnBdry(e))
        PCU_ALWAYS_ASSERT(data[i] == pumi_ment_getGlobalID(e)*(1+pumi_ment_getNumRmt(e)));
      else
        PCU_ALWAYS_ASSERT(data[i] == xyz[i]);
  }
  m->end(it);

  pumi_ghost_delete(m);

  for (int i=0; i<4; ++i)
  {
    if (org_mcount[i] != pumi_mesh_getNumEnt(m, i))
       std::cout<<"("<<pumi_rank()<<") ERROR dim "<<i<<": org ent count "<<org_mcount[i]<<", current ent count "<<pumi_mesh_getNumEnt(m, i)<<"\n";
    PCU_ALWAYS_ASSERT(org_mcount[i] == pumi_mesh_getNumEnt(m, i));
  }
  
  delete [] org_mcount;
}
