/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <mpi.h>
#include <parma.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfZoltan.h>
#include <pcu_util.h>
#include <iostream>
#include <string.h>
#include <map>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>
#include "apf.h"
#include "apfShape.h"
#include "apfNumbering.h"

using std::map;

// mesh creation
pMesh pumi_mesh_create(pGeom g, int mesh_dim, bool periodic)
{
  pumi::instance()->mesh = apf::makeEmptyMdsMesh(g->getGmi(), mesh_dim, periodic);
  return pumi::instance()->mesh;
}

void pumi_mesh_freeze(pMesh m)
{
  deriveMdsModel(m);
  m->acceptChanges();
  pumi_geom_freeze(pumi_mesh_getGeom(m));
}

pMeshEnt pumi_mesh_createVtx(pMesh m, pGeomEnt ge, double* xyz)
{
  apf::Vector3 coord(xyz[0],xyz[1],xyz[2]);
  apf::Vector3 param(0,0,0);
  return m->createVertex((apf::ModelEntity*)ge, coord, param);
}

pMeshEnt pumi_mesh_createEnt(pMesh m, pGeomEnt ge, int ent_topology, pMeshEnt* down)
{
  return m->createEntity(ent_topology, (apf::ModelEntity*)ge, down);
}

pMeshEnt pumi_mesh_createElem(pMesh m, pGeomEnt ge, int ent_topology, pMeshEnt* vertices) 
{
  return apf::buildElement(m, (apf::ModelEntity*)ge, ent_topology, vertices);
}

void generate_globalid(pMesh m, pMeshTag tag, int dim, pOwnership o)
{
  pMeshEnt e;
  int num_own=0;

  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    if ((o && o->isOwned(e)) || (!o && m->isOwned(e)))
      ++num_own;
  m->end(it);

  PCU_Exscan_Ints(&num_own,1);
  int initial_id=num_own;

  PCU_Comm_Begin();
  it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    if ((o && !o->isOwned(e)) || (!o && !m->isOwned(e)))
      continue;

    m->setIntTag(e, tag, &initial_id);
    Copies remotes;
    m->getRemotes(e, remotes);
    APF_ITERATE(Copies, remotes, it)
    {
      PCU_COMM_PACK(it->first, it->second);
      PCU_Comm_Pack(it->first, &initial_id, sizeof(int));
    }

    if (m->isGhosted(e))
    {
      Copies ghosts;
      m->getGhosts(e, ghosts);
      APF_ITERATE(Copies, ghosts, it)
      {
        PCU_COMM_PACK(it->first, it->second);
        PCU_Comm_Pack(it->first, &initial_id, sizeof(int));
      }
    }
    ++initial_id;
  }
  m->end(it);

  PCU_Comm_Send();
  int global_id;
  while (PCU_Comm_Listen())
    while (!PCU_Comm_Unpacked())
    {
      pMeshEnt remote_ent;
      PCU_COMM_UNPACK(remote_ent);
      PCU_Comm_Unpack(&global_id, sizeof(int));
      m->setIntTag(remote_ent, tag, &global_id);
    }
}

//*******************************************************
void pumi_mesh_createGlobalID(pMesh m, pOwnership o)
//*******************************************************
{
  pMeshTag tag = m->findTag("global_id");
  if (tag)  // destroy existing tag
  {
    for (int i=0; i<4; ++i)
      apf::removeTagFromDimension(m, tag, m->getDimension());
  }  
  else
    tag = m->createIntTag("global_id",1);

  for (int i=0; i<=m->getDimension(); ++i)
    generate_globalid(m, tag, i, o);
}

//*******************************************************
void pumi_mesh_deleteGlobalID(pMesh m)
//*******************************************************
{
  pMeshTag tag = m->findTag("global_id");
  if (!tag) return;

  for (int i=0; i<4; ++i)
    apf::removeTagFromDimension(m, tag, i);

  m->destroyTag(tag);
}

// *********************************************************
pumi::pumi(): mesh(NULL), model(NULL) 
{
  ghost_tag=NULL;
  ghosted_tag=NULL;
  num_local_ent = NULL;
  num_own_ent = NULL;
  num_global_ent = NULL;
}

pumi::~pumi()
{
  if (num_own_ent) 
  {
    delete [] num_local_ent;
    delete [] num_own_ent;
    delete [] num_global_ent;
  }
}



apf::Migration* getPlan(apf::Mesh* m, int num_target_part)
{
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  pMeshTag weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, num_target_part);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void split_comm(int num_out_comm)
{
  int self = PCU_Comm_Self();
  int group_id = self % num_out_comm;
  int in_group_rank = self / num_out_comm;
  MPI_Comm groupComm;
  MPI_Comm_split(PCU_Get_Comm(), group_id, in_group_rank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void merge_comm(MPI_Comm oldComm)
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(oldComm);
  MPI_Comm_free(&prevComm);
}


pGeom pumi_mesh_getGeom(pMesh)
{
  return pumi::instance()->model;
}

// load a serial mesh on master process then distribute as per the distribution object
pMesh pumi_mesh_loadSerial(pGeom g, const char* filename, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  MPI_Comm prevComm = PCU_Get_Comm();
  int num_target_part = PCU_Comm_Peers();
  bool isMaster = ((PCU_Comm_Self() % num_target_part) == 0);
  pMesh m = 0;
  split_comm(num_target_part);
  if (isMaster) 
    m = apf::loadMdsMesh(g->getGmi(), filename);
  merge_comm(prevComm);
  pumi::instance()->mesh = expandMdsMesh(m, g->getGmi(), 1);
  return pumi::instance()->mesh;
}

pMesh pumi_mesh_load(pMesh m)
{ //TODO add sanity checks
  pumi::instance()->mesh = m;
  pumi_mesh_print(pumi::instance()->mesh);
  return pumi::instance()->mesh;
}


pMesh pumi_mesh_load(pGeom g, const char* filename, int num_in_part, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  if (num_in_part==1 && pumi_size()>1) // do static partitioning
  {
    MPI_Comm prevComm = PCU_Get_Comm();
    int num_target_part = PCU_Comm_Peers()/num_in_part;
    bool isMaster = ((PCU_Comm_Self() % num_target_part) == 0);
    pMesh m = 0;
    apf::Migration* plan = 0;   
    split_comm(num_target_part);
    if (isMaster) {
      m = apf::loadMdsMesh(g->getGmi(), filename);
      plan = getPlan(m, num_target_part);
    }
    merge_comm(prevComm);
    pumi::instance()->mesh = apf::repeatMdsMesh(m, g->getGmi(), plan, num_target_part);
  }
  else
    pumi::instance()->mesh = apf::loadMdsMesh(g->getGmi(), filename);
  pumi_mesh_print(pumi::instance()->mesh);
  return pumi::instance()->mesh;
}



void send_entities(pMesh mesh, int dim)
{
  int local_id, self = PCU_Comm_Self();
  pMeshEnt e;
  pMeshIter it = mesh->begin(dim);
  while ((e = mesh->iterate(it)))
  {
    local_id = getMdsIndex(mesh, e);
    for (int pid=0; pid<PCU_Comm_Peers(); ++pid)
    {
      if (pid==self) continue;
      PCU_Comm_Pack(pid, &local_id, sizeof(int));
      PCU_COMM_PACK(pid, e);
    }
  }
  mesh->end(it);
}

#include "apfMDS.h"
#include "apfPM.h"
pMesh pumi_mesh_loadAll(pGeom g, const char* filename, bool stitch_link)
{
  if (pumi_size()==1) 
    pumi::instance()->mesh = apf::loadMdsMesh(g->getGmi(), filename);
  else
  {
    double t0 = PCU_Time();
    MPI_Comm prevComm = PCU_Get_Comm();
    int num_target_part = PCU_Comm_Peers();
    split_comm(num_target_part);
    // no pmodel & remote links setup
    pumi::instance()->mesh = apf::loadSerialMdsMesh(g->getGmi(), filename); 
    merge_comm(prevComm);
    if (!PCU_Comm_Self())
      lion_oprint(1,"serial mesh %s loaded in %f seconds\n", filename, PCU_Time() - t0);
  }

  if (pumi_size()>1 && stitch_link) 
  {
    stitchMesh(pumi::instance()->mesh);
    pumi::instance()->mesh->acceptChanges();
  }

  return pumi::instance()->mesh;
}

void pumi_mesh_migrate(pMesh m, Migration* plan)
{
  apf::migrate(m, plan);
}

int pumi_mesh_getDim(pMesh m)
{
  return m->getDimension();
}

void pumi_mesh_setCount(pMesh m, pOwnership o)
{
#ifdef DEBUG
  if (!pumi_rank()) std::cout<<"[PUMI INFO] "<<__func__<<" start\n";
#endif

  if (!pumi::instance()->num_local_ent)
  { 
    pumi::instance()->num_local_ent = new int[4];
    pumi::instance()->num_own_ent = new int[4];
    pumi::instance()->num_global_ent = new int[4];
  }

  for (int dim=0; dim<4; ++dim)
  {
    pumi::instance()->num_local_ent[dim]=m->count(dim);
    if (!o) // NULL
      pumi::instance()->num_own_ent[dim] = countOwned(m, dim);
    else
    {
      apf::MeshIterator* it = m->begin(dim);
      apf::MeshEntity* e;
      int n = 0;
      while ((e = m->iterate(it)))
        if (o->isOwned(e))
          ++n;
      m->end(it);
      pumi::instance()->num_own_ent[dim] = n;
    }
  }
  MPI_Allreduce(pumi::instance()->num_own_ent, pumi::instance()->num_global_ent, 4, MPI_INT, MPI_SUM, PCU_Get_Comm());
#ifdef DEBUG
  if (!pumi_rank()) std::cout<<"[PUMI INFO] "<<__func__<<" end\n";
#endif
}

int pumi_mesh_getNumEnt(pMesh m, int dim)
{ return m->count(dim); }

int pumi_mesh_getNumOwnEnt(pMesh m, int dim)
{ 
  PCU_ALWAYS_ASSERT(pumi::instance()->num_own_ent);
  if (pumi::instance()->num_local_ent[dim]!=(int)m->count(dim) && !PCU_Comm_Self()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<": mesh count is not set. Please call pumi_mesh_setCount\n";
    return -1;
  }
  return pumi::instance()->num_own_ent[dim];
}

int pumi_mesh_getNumGlobalEnt(pMesh m, int dim)
{ 
  PCU_ALWAYS_ASSERT(pumi::instance()->num_global_ent);
  if (pumi::instance()->num_local_ent[dim]!=(int)m->count(dim) && !PCU_Comm_Self()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<": mesh count is not set. Please call pumi_mesh_setCount\n";
    return -1;
  }
  return pumi::instance()->num_global_ent[dim];
}

pMeshEnt pumi_mesh_findEnt(pMesh m, int d, int id)
{
  return getMdsEntity(m, d, id);
}

#include <parma.h>

void print_copies(pMesh m, pMeshEnt e)
{
  if (m->isShared(e))
  {
    Copies remotes;
    m->getRemotes(e,remotes);
    std::cout<<"\tremotes: ";
    APF_ITERATE(Copies,remotes,rit)
      std::cout<<"("<<rit->first<<", "<<rit->second<<") ";
    std::cout<<"\n";
  }
  if (m->isGhosted(e) || m->isGhost(e))
  {
    Copies ghosts;
    m->getGhosts(e,ghosts);
    std::cout<<"\tghosts: ";
    APF_ITERATE(Copies,ghosts,rit)
      std::cout<<"("<<rit->first<<", "<<rit->second<<") ";
    std::cout<<"\n";
  }
}

void pumi_mesh_print (pMesh m, bool print_ent)
{
  if (!PCU_Comm_Self()) std::cout<<"\n=== mesh size and tag info === \n";

  int* local_entity_count = new int[4*PCU_Comm_Peers()];
  int* own_entity_count = new int[4*PCU_Comm_Peers()];

  for (int i=0; i<4*PCU_Comm_Peers();++i)
    local_entity_count[i]=own_entity_count[i]=0;

  pMeshEnt e;
  int self = pumi_rank();

  for (int d=0; d<4;++d)
  {
    local_entity_count[4*self+d] = m->count(d);
    pMeshIter it = m->begin(d); // face
    while ((e = m->iterate(it)))
    {
      if (m->getOwner(e)==self)
        ++own_entity_count[4*pumi_rank()+d];
    }
    m->end(it);
  }
  
  int* global_local_entity_count = new int[4*PCU_Comm_Peers()]; 
  int* global_own_entity_count = new int[4*PCU_Comm_Peers()]; 

  MPI_Allreduce(local_entity_count, global_local_entity_count, 4*PCU_Comm_Peers(), 
                MPI_INT, MPI_SUM, PCU_Get_Comm());

  MPI_Allreduce(own_entity_count, global_own_entity_count, 4*PCU_Comm_Peers(), 
                MPI_INT, MPI_SUM, PCU_Get_Comm());

 
  if (!PCU_Comm_Self())
  {
    int* global_entity_count = new int[4]; 
    global_entity_count[0]=global_entity_count[1]=global_entity_count[2]=global_entity_count[3]=0;
    for (int d=0; d<4;++d)
    {
      for (int p=0; p<PCU_Comm_Peers();++p)
        global_entity_count[d] += global_own_entity_count[p*4+d];
    }

     std::cout<<"# global ent: v "<<global_entity_count[0]<<", e "<<global_entity_count[1]
              <<", f "<<global_entity_count[2]<<", r "<<global_entity_count[3]<<"\n";

    delete [] global_entity_count;

    for (int p=0; p<PCU_Comm_Peers(); ++p)
      std::cout<<"(p"<<p<<") # local ent: v "<<global_local_entity_count[p*4]
        <<", e "<<global_local_entity_count[p*4+1]
        <<", f "<<global_local_entity_count[p*4+2]
        <<", r "<<global_local_entity_count[p*4+3]<<"\n";
    for (int p=0; p<PCU_Comm_Peers(); ++p)
      if (global_own_entity_count[p*4])
        std::cout<<"(p"<<p<<") # own ent: v "<<global_own_entity_count[p*4]
          <<", e "<<global_own_entity_count[p*4+1]
          <<", f "<<global_own_entity_count[p*4+2]
          <<", r "<<global_own_entity_count[p*4+3]<<"\n";
    std::cout<<"\n";
  }

  delete [] local_entity_count;
  delete [] global_local_entity_count;
  delete [] own_entity_count;
  delete [] global_own_entity_count;

  if (!PCU_Comm_Self()) 
  {
    std::cout<<"mesh shape: \""<< m->getShape()->getName()<<"\"\n";

    apf::DynamicArray<pMeshTag> tags; // tags
    m->getTags(tags);
    int n = tags.getSize();
    if (n) 
      for (int i = 0; i < n; ++i) 
        std::cout<<"tag "<<i<<": \""<< m->getTagName(tags[i])<<"\", type "
                << m->getTagType(tags[i])<<", size "<< m->getTagSize(tags[i])<<"\n";

    if (m->countFields())
      for (int i = 0; i < m->countFields(); ++i)  // fields
        std::cout<<"field "<<i<<": \""<< getName(m->getField(i))
                 <<"\", size "<<  apf::countComponents(m->getField(i))<<"\n";

    if (m->countNumberings()) 
      for (int i = 0; i < m->countNumberings(); ++i)  // fields
        std::cout<<"numbering "<<i<<": \""<< getName(m->getNumbering(i))
                 <<"\" on shape \""<< getShape(m->getNumbering(i))->getName()<<"\n";
  
    if (m->countGlobalNumberings()) 
      for (int i = 0; i < m->countGlobalNumberings(); ++i)  // fields
        std::cout<<"numbering "<<i<<": \""<< getName(m->getGlobalNumbering(i))
                 <<"\" on shape \""<< getShape(m->getGlobalNumbering(i))->getName()<<"\n";
  }

  // print mesh entities
  if (!print_ent) return;

  if (!m->findTag("global_id")) 
  {  
    pumi_mesh_createGlobalID(m);
    if (!PCU_Comm_Self()) std::cout<<__func__<<": global id generated\n";
  }

  int global_id;
  apf::MeshIterator* vit = m->begin(0);
  while ((e = m->iterate(vit)))
  {
    apf::Vector3 xyz;
    m->getPoint(e, 0, xyz);
    if (m->isGhost(e))
      std::cout<<"("<<PCU_Comm_Self()<<") GHOST vtx "<<e<<": "<<pumi_ment_getGlobalID(e)
             <<" ("<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<")\n";
    else
      std::cout<<"("<<PCU_Comm_Self()<<") vtx "<<e<<": "<<pumi_ment_getGlobalID(e)
             <<" ("<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<")\n";
    print_copies(m,e);
  }
  m->end(vit);

  for (int d=1; d<4; ++d)
  {
    apf::MeshIterator* eit = m->begin(d);
    while ((e = m->iterate(eit)))
    {
      global_id=pumi_ment_getGlobalID(e);
      apf::Downward down;
      int num_down=m->getDownward(e,d-1,down); 
      if (m->isGhost(e)) 
        std::cout<<"("<<PCU_Comm_Self()<<") GHOST e "<<e<<": [d "<<d<<", id "<<global_id<<"] down: ";
      else
        std::cout<<"("<<PCU_Comm_Self()<<") e "<<e<<": [d "<<d<<", id "<<global_id<<"] down: ";
      for (int i=0; i<num_down; ++i)
        std::cout<<pumi_ment_getGlobalID(down[i])<<" ";
      std::cout<<"\n";
      print_copies(m,e);
    }
    m->end(eit);
  }
  pumi_mesh_deleteGlobalID(m);
}

void pumi_mesh_write (pMesh m, const char* filename, const char* mesh_type)
{
  if (!strcmp(mesh_type,"mds"))
    m->writeNative(filename);
  else if (!strcmp(mesh_type,"vtk"))
  {
    // attach ghost flag for visualization
    apf::Field* ghost_f = apf::createStepField(m, "ghost_field", apf::SCALAR);
    apf::Field* own_f = apf::createStepField(m, "own_field", apf::SCALAR);
    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    double ghost_value, own_partid;
    int ghost_flag=1, non_ghost_flag=0;
    while ((e = m->iterate(it))) 
    { 
      own_partid=pumi_ment_getOwnPID(e);
      if (m->isGhost(e))
        ghost_value=ghost_flag;
      else
        ghost_value=non_ghost_flag;
      setScalar(ghost_f, e, 0, ghost_value);
      setScalar(own_f, e, 0, own_partid);
    }
    m->end(it);

    apf::writeVtkFiles(filename, m);

    destroyField(ghost_f);
    destroyField(own_f);
  }
  else
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
}

void pumi_mesh_delete(pMesh m)
{
  if (m->findTag("ghost_tag"))
    m->destroyTag(pumi::instance()->ghost_tag);
  if (m->findTag("ghosted_tag"))
    m->destroyTag(pumi::instance()->ghosted_tag);
  m->destroyNative();
  apf::destroyMesh(m);
}

// create/delete direct Adjacency for all entities except for one-level apart
bool pumi_mesh_hasAdjacency(pMesh m, int from_dim, int to_dim)
{
  return m->hasAdjacency(from_dim, to_dim);
}

void pumi_mesh_createAdjacency(pMesh m, int from_dim, int to_dim)
{
  m->createAdjacency(from_dim, to_dim);
}

void pumi_mesh_deleteAdjacency(pMesh m, int from_dim, int to_dim)
{
  m->deleteAdjacency(from_dim, to_dim);
}

void pumi_mesh_createFullAdjacency(pMesh m)
{
  if (m->getDimension()==3)
  {
    pumi_mesh_createAdjacency(m,3,1);
    pumi_mesh_createAdjacency(m,1,3);
    pumi_mesh_createAdjacency(m,3,0);
    pumi_mesh_createAdjacency(m,0,3);
  }
  pumi_mesh_createAdjacency(m,2,0);
  pumi_mesh_createAdjacency(m,0,2);
}

void pumi_mesh_verify(pMesh m, bool abort_on_error)
{
  apf::verify(m, abort_on_error);
}

#include <assert.h>

void pumi_ownership_verify(pMesh m, pOwnership o)
{
  if (!pumi::instance()->num_local_ent) {
    pumi_mesh_setCount(m, o);
  }

    int num_own_ent, own_partid, mesh_dim = m->getDimension()-1;
    pMeshEnt e;
    pMeshEnt own_copy;
    pMeshEnt remote_copy;

    for (int d=0; d<mesh_dim+1; ++d)
    {
      num_own_ent=0;
      pMeshIter it = m->begin(d);
      while ((e = m->iterate(it)))
      {
          own_partid = pumi_ment_getOwnPID(e,o);
          own_copy = pumi_ment_getOwnEnt(e,o);
          if (!own_copy)
          {
            std::cout<<"[ERROR] ("<<pumi_rank()<<") "<<__func__<<": pumi_ment_getOwnEnt(dim "<<d<<", id "<<pumi_ment_getID(e)<<", pid "<<own_partid<<") not found\n";
            print_copies(m,e);
          }
          assert(own_copy);
          if (own_partid==PCU_Comm_Self()) 
          {
            ++num_own_ent;
            if (own_copy!=e)
            {
              std::cout<<"[ERROR] ("<<pumi_rank()<<") "<<__func__<<": pumi_ment_getOwnEnt(dim "<<d<<", id "<<pumi_ment_getID(e)<<") is not self copy\n";
              print_copies(m,e);
            }
            assert(own_copy==e);
          }
          else
          {
            remote_copy = pumi_ment_getRmt(e, own_partid); 
            if (!remote_copy)
            {
              std::cout<<"[ERROR] ("<<pumi_rank()<<") "<<__func__<<": pumi_ment_getRmt(dim "<<d<<", id "<<pumi_ment_getID(e)<<", pid "<<own_partid<<") not found\n";
              print_copies(m,e);
            }
            if (own_copy!=remote_copy)
            {
              std::cout<<"[ERROR] ("<<pumi_rank()<<") "<<__func__<<": pumi_ment_getRmt and pumi_ment_getOwnEnt mismatch for e(dim "<<d<<", id "<<pumi_ment_getID(e)<<", pid "<<own_partid<<")\n";
              print_copies(m,e);
            }
            assert(own_copy==remote_copy);
          }
      }
      m->end(it);
      if (pumi_mesh_getNumOwnEnt(m, d) != num_own_ent)
         std::cout<<"[ERROR] "<<__func__<<": pumi_mesh_getNumOwnEnt(m, "<<d<<") "<<pumi_mesh_getNumOwnEnt(m, d)<<", num_own_ent "<<num_own_ent<<" -- check isOwned!\n";
      assert (pumi_mesh_getNumOwnEnt(m, d) == num_own_ent);
    } // for

  if (!pumi_rank()) std::cout<<__func__<<": ownership is valid\n";
}

Distribution::Distribution(pMesh mesh)
{
  m = mesh;
  parts_vec = NULL;
  element_count=0;
}

Distribution::~Distribution()
{
  if (parts_vec)
    delete [] parts_vec;
}


bool Distribution::has(pMeshEnt e)
{
  int i = getMdsIndex(m, e);
  if (parts_vec[i].size())
    return true;
  else
    return false;
}

void Distribution::send(pMeshEnt e, int to)
{
  if (parts_vec == NULL) 
  {
    int dim = m->getDimension();
    int nele = m->count(dim);
    parts_vec = new Parts[nele];
  }
  int i = getMdsIndex(m, e);
  parts_vec[i].insert(to);
}

Parts& Distribution::sending(pMeshEnt e)
{
  int i = getMdsIndex(m, e);
  PCU_ALWAYS_ASSERT(parts_vec[i].size()>0);
  return parts_vec[i];
}

int Distribution::count()
{
  if (element_count==0)
  {
    int nele = m->count(m->getDimension());
    for (int i=0; i<nele; ++i)
      if (parts_vec[i].size()>0)
        ++element_count;
  }
  return element_count;
}

void Distribution::print()
{
  pMeshEnt e;
  apf::MeshIterator* it = m->begin(m->getDimension());
  int i=-1;
  while ((e = m->iterate(it)))
  {
    ++i;
    if (parts_vec[i].size()==0) continue;
    APF_ITERATE(Parts,parts_vec[i],pit)
      std::cout<<"("<<PCU_Comm_Self()<<") distribute element "<<i<<" to "<<*pit<<"\n";

  }
  m->end(it);
}

static void distr_getAffected (pMesh m, Distribution* plan, EntityVector affected[4])
{
  int maxDimension = m->getDimension();
  affected[maxDimension].reserve(plan->count());

  pMeshEnt e;
  apf::MeshIterator* it = m->begin(maxDimension);
  int i=0;
  while ((e = m->iterate(it)))
  {
    if (plan->parts_vec[i].size()>0) 
      affected[maxDimension].push_back(e);
    ++i;
  }
  m->end(it);

  int dummy=1;
  pMeshTag tag = m->createIntTag("distribution_affected",1);
  for (int dimension=maxDimension-1; dimension >= 0; --dimension)
  {
    int upDimension = dimension + 1;
    PCU_Comm_Begin();
    APF_ITERATE(EntityVector,affected[upDimension],it)
    {
      pMeshEnt up = *it;
      apf::Downward adjacent;
      int na = m->getDownward(up,dimension,adjacent);
      for (int i=0; i < na; ++i)
      {
        if ( ! m->hasTag(adjacent[i],tag))
        {
          m->setIntTag(adjacent[i],tag,&dummy);
          affected[dimension].push_back(adjacent[i]);
        }
        Copies remotes;
        m->getRemotes(adjacent[i],remotes);
        APF_ITERATE(Copies,remotes,rit)
          PCU_COMM_PACK(rit->first,rit->second);
        if (m->hasMatching())
        {
          apf::Matches matches;
          m->getMatches(adjacent[i],matches);
          for (size_t j=0; j < matches.getSize(); ++j)
            PCU_COMM_PACK(matches[j].peer,matches[j].entity);
        }
      }//downward adjacent loop
    }//upward affected loop
    PCU_Comm_Send();
    while (PCU_Comm_Receive())
    {
      pMeshEnt entity;
      PCU_COMM_UNPACK(entity);
      if ( !m->hasTag(entity,tag))
      {
        m->setIntTag(entity,tag,&dummy);
        affected[dimension].push_back(entity);
      }
    }
    APF_ITERATE(EntityVector,affected[dimension],it)
      m->removeTag(*it,tag);
  }//dimension loop
  m->destroyTag(tag);
}

static Parts distr_makeResidence(Parts& parts)
{
  Parts r;
  APF_ITERATE(Parts, parts, pit)
    r.insert(*pit);
  return r;
}

// *********************************************************
static void distr_updateResidences(pMesh m,
    Distribution* plan, EntityVector affected[4])
// *********************************************************
{
  int maxDimension = m->getDimension();
  pMeshEnt e;

  apf::MeshIterator* it = m->begin(maxDimension);
  int i=0;
  while ((e = m->iterate(it)))
  {
    if (plan->parts_vec[i].size()>0) 
    {
      Parts res = distr_makeResidence(plan->parts_vec[i]);
      m->setResidence(e,res);
    }
    ++i;
  }
  m->end(it);

  for (int dimension = maxDimension-1; dimension >= 0; --dimension)
  {
    PCU_Comm_Begin();
    APF_ITERATE(EntityVector,affected[dimension],it)
    {
      pMeshEnt entity = *it;
      Parts newResidence;
      apf::Up upward;
      m->getUp(entity, upward);
      for (int ui=0; ui < upward.n; ++ui)
      {
        pMeshEnt up = upward.e[ui];
        Parts upResidence;
        m->getResidence(up,upResidence);
        apf::unite(newResidence,upResidence);
      }
      m->setResidence(entity,newResidence);
      Copies remotes;
      m->getRemotes(entity,remotes);
      APF_ITERATE(Copies,remotes,rit)
      {
        PCU_COMM_PACK(rit->first,rit->second);
        apf::packParts(rit->first,newResidence);
      }
    }
    PCU_Comm_Send();
    while(PCU_Comm_Receive())
    {
      pMeshEnt entity;
      PCU_COMM_UNPACK(entity);
      Parts current;
      m->getResidence(entity,current);
      Parts incoming;
      apf::unpackParts(incoming);
      apf::unite(current,incoming);
      m->setResidence(entity,current);
    }
  }
}

// *********************************************************
void distribute(pMesh m, Distribution* plan)
// *********************************************************
{
  EntityVector affected[4];
  distr_getAffected(m,plan,affected);
  EntityVector senders[4];
  getSenders(m,affected,senders);
  reduceMatchingToSenders(m,senders);
  distr_updateResidences(m,plan,affected);
  delete plan;
  moveEntities(m,senders);
  updateMatching(m,affected,senders);
  deleteOldEntities(m,affected);
  m->acceptChanges();
}

// *********************************************************
void pumi_mesh_distribute(pMesh m, Distribution* plan)
// *********************************************************
{
  if (PCU_Comm_Peers()==1) return;

  if (pumi::instance()->ghosted_tag)
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" not supported with ghosted mesh\n";
    return;
  }
  distribute(m, plan);
}
