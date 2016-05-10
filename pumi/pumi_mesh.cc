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

apf::Migration* getPlan(apf::Mesh* m, int partitionFactor)
{
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void switchToOriginals(int partitionFactor)
{
  int self = PCU_Comm_Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}

pMesh pumi_mesh_create(pGeom g, const char* filename, int num_in_part, int num_proc_grp, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }

  if (num_in_part==PCU_Comm_Peers())
    pumi::instance()->mesh = apf::loadMdsMesh(g, filename);
  else
  {
  // if do static partitioning
  int partitionFactor = PCU_Comm_Peers()/num_in_part;
  bool isOriginal = ((PCU_Comm_Self() % partitionFactor) == 0);
  pMesh m = 0;
  apf::Migration* plan = 0;
  switchToOriginals(partitionFactor);
  if (isOriginal) {
    m = apf::loadMdsMesh(g, filename);
    plan = getPlan(m, partitionFactor);
  }
  switchToAll();
  pumi::instance()->mesh = apf::repeatMdsMesh(m, g, plan, partitionFactor);
  }
  return pumi::instance()->mesh;
}

int pumi_mesh_getdim(pMesh m)
{
  return m->getDimension();
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


