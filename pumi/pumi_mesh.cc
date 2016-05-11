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

// *********************************************************
pumi::pumi()
// *********************************************************
{
  mesh = NULL;
  model = NULL;
}

// *********************************************int num_proc_grp************
pumi::~pumi()
// *********************************************************
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

  return pumi::instance()->mesh;
}

int pumi_mesh_getdim(pMesh m)
{
  return m->getDimension();
}

#include <parma.h>
void pumi_mesh_print (pMesh m)
{
  if (!PCU_Comm_Self()) std::cout<<"=== mesh info === \n";
  printStats(m);
  Parma_PrintPtnStats(m, "initial");
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
  m->destroyNative();
  apf::destroyMesh(m);
}


