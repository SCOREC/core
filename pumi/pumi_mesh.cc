/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
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
#include <assert.h>
#include <iostream>
#include <string.h>

//*******************************************************
void destroy_node_global_numbering(apf::Mesh2* m)
//*******************************************************
{
  if (m->findField("node own partid field"))
    destroyField(m->findField("node own partid field"));
  if (m->findField("node global id field"))
    destroyField(m->findField("node global id field"));
}

//*******************************************************
void generate_node_global_numbering(apf::Mesh2* m)
//*******************************************************
{
//  if (!PCU_Comm_Self())  std::cout<<"[M3D-C1 INFO] ***** GENERATING GLOBAL NODE NUMBERING ***** \n"; 
  destroy_node_global_numbering(m);

  double id[1];
  double own_partid[1];
  apf::Field* node_ownpid_f = createPackedField(m, "node own partid field", 1);
  apf::freeze(node_ownpid_f);
  apf::Field* node_globalid_f = createPackedField(m, "node global id field", 1);
  apf::freeze(node_globalid_f);

  // count #own_vtx
  int num_own_ent=0;
  pMeshEnt e;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (m->getOwner(e)==PCU_Comm_Self())
      ++num_own_ent;
  }
  m->end(it);

  // generate global node_id
  pumi::instance()->num_own_vtx=num_own_ent;
  PCU_Exscan_Ints(&num_own_ent,1);
  int start=num_own_ent;

  PCU_Comm_Begin();

  it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid[0]=(double)m->getOwner(e); 
    setComponents(node_ownpid_f, e, 0, own_partid);    
    if ((int)(own_partid[0])!=PCU_Comm_Self()) continue;
    id[0] = (double) start;
    setComponents(node_globalid_f, e, 0, id);
    apf::Copies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(apf::Copies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&start,sizeof(int));
    }
    ++start;
  }
  m->end(it);
  PCU_Comm_Send();

  int value;
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&value,sizeof(int));
      id[0] = (double) value;
      setComponents(node_globalid_f, r, 0, id);
    }
}


// *********************************************************
pumi::pumi(): mesh(NULL), model(NULL), org_mesh(NULL) {}

pumi::~pumi()
{
  delete _instance;
  _instance = NULL;
}

pumi* pumi::_instance=NULL;
pumi* pumi::instance()
{
  if (_instance==NULL)
    _instance = new pumi();
  return _instance;
}

apf::Migration* getPlan(apf::Mesh* m, int num_target_part)
{
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
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

void merge_comm()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}

pMesh pumi_mesh_create(pGeom g, const char* filename, int num_in_part, 
                       int num_proc_group, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  if (num_proc_group==1)
  {
    if (num_in_part==PCU_Comm_Peers())
      pumi::instance()->mesh = apf::loadMdsMesh(g, filename);
    else if (num_in_part==1) // do static partitioning
    {
      int num_target_part = PCU_Comm_Peers()/num_in_part;
      bool isMaster = ((PCU_Comm_Self() % num_target_part) == 0);
      pMesh m = 0;
      apf::Migration* plan = 0;
      split_comm(num_target_part);
      if (isMaster) {
        m = apf::loadMdsMesh(g, filename);
        plan = getPlan(m, num_target_part);
      }
      merge_comm();
      pumi::instance()->mesh = apf::repeatMdsMesh(m, g, plan, num_target_part);
    }
  }
  else // load mesh in each process group 
  {
    assert(PCU_Comm_Peers()/num_in_part==num_proc_group);
    //split_comm(PCU_Comm_Peers()/num_proc_group);
  
    int self = PCU_Comm_Self();
    int group_size=PCU_Comm_Peers()/num_proc_group;
    int group_id = self/group_size; // divide
    int in_group_rank = self%group_size; // modulo
    MPI_Comm groupComm;
    MPI_Comm_split(PCU_Get_Comm(), group_id, in_group_rank, &groupComm);
    PCU_Switch_Comm(groupComm);

    // each proc group loads mesh  
    pumi::instance()->mesh = apf::loadMdsMesh(g, filename);
//    apf::disownMdsModel(m3dc1_mesh::instance()->mesh);
    merge_comm();
    // FIXME: in non-master process group, incorrect proc ranks in 
    //        remote copies, owning part and partition model entity's residence part 
  }
  generate_node_global_numbering(pumi::instance()->mesh);
  return pumi::instance()->mesh;
}

int pumi_mesh_getdim(pMesh m)
{
  return m->getDimension();
}

int pumi_mesh_getnument(pMesh m, int dim)
{ return m->count(dim); }

#include <parma.h>
void pumi_mesh_print (pMesh m)
{
  if (!PCU_Comm_Self()) std::cout<<"\n=== mesh count === \n";
  printStats(m);
  for (int i=0; i<PCU_Comm_Peers(); ++i)
  {
    if (i==pumi_rank())
      std::cout<<"(p"<<PCU_Comm_Self()<<") # local ent: v "<<m->count(0)
        <<", e "<<m->count(1)<<", f "<<m->count(2)<<", r "<<m->count(3)<<"\n";
    MPI_Barrier(MPI_COMM_WORLD);
  }
  //Parma_PrintPtnStats(m, "initial");
}

void pumi_mesh_write (pMesh m, const char* filename, const char* mesh_type)
{
  if (!strcmp(mesh_type,"mds"))
    m->writeNative(filename);
  else if (!strcmp(mesh_type,"vtk"))
    apf::writeVtkFiles(filename, m);
  else
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
}

void pumi_mesh_delete(pMesh m)
{
  destroy_node_global_numbering(pumi::instance()->mesh);
  pumi::instance()->mesh->destroyNative();
  apf::destroyMesh(pumi::instance()->mesh);
}


