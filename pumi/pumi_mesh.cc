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
#include <assert.h>
#include <cstdlib>
#include "apf.h"

//*******************************************************
void destroy_global_numbering(apf::Mesh2* m)
//*******************************************************
{
  pMeshTag tag = m->findTag("global_id");
  if (!tag) return;
  for (int i=0; i<4; ++i)
    apf::removeTagFromDimension(m, tag, m->getDimension());

  m->destroyTag(tag);
}

// *********************************************************
pumi::pumi(): mesh(NULL), model(NULL) 
{
  ghost_tag=NULL;
  ghosted_tag=NULL;
}

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


pMesh pumi_mesh_load(pGeom g, const char* filename, int num_in_part, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  if (num_in_part==1) // do static partitioning
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

int pumi_mesh_getNumEnt(pMesh m, int dim)
{ return m->count(dim); }

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

void pumi_mesh_print (pMesh m, int p)
{
  if (!PCU_Comm_Self()) std::cout<<"\n=== mesh size and tag info === \nglobal ";
  printStats(m);

  int* local_entity_count = new int[4*PCU_Comm_Peers()];
  for (int i=0; i<4*PCU_Comm_Peers();++i)
    local_entity_count[i]=0;
  for (int d=0; d<4;++d)
    local_entity_count[4*pumi_rank()+d] = m->count(d);
  
  int* global_entity_count = new int[4*PCU_Comm_Peers()]; 

  MPI_Allreduce(local_entity_count, global_entity_count, 4*PCU_Comm_Peers(), MPI_INT, MPI_SUM, PCU_Get_Comm());
  if (!PCU_Comm_Self())
  {
    for (int p=0; p<PCU_Comm_Peers(); ++p)
      std::cout<<"(p"<<p<<") # local ent: v "<<global_entity_count[p*4]
        <<", e "<<global_entity_count[p*4+1]
        <<", f "<<global_entity_count[p*4+2]
        <<", r "<<global_entity_count[p*4+3]<<"\n";
  }

  delete [] local_entity_count;
  delete [] global_entity_count;

  if (!PCU_Comm_Self())
  {
    apf::DynamicArray<pMeshTag> tags;
    m->getTags(tags);
    int n = tags.getSize();
    for (int i = 0; i < n; ++i) 
      std::cout<<"tag "<<i<<": \""<< m->getTagName(tags[i])<<"\", type "
              << m->getTagType(tags[i])<<", size "<< m->getTagSize(tags[i])<<"\n";
  }

  // print mesh entities
  if (p!=PCU_Comm_Self()) return;

  pMeshEnt e;
  int global_id;
  apf::MeshIterator* vit = m->begin(0);
  while ((e = m->iterate(vit)))
  {
    apf::Vector3 xyz;
    m->getPoint(e, 0, xyz);
    if (m->isGhost(e))
      std::cout<<"("<<PCU_Comm_Self()<<") GHOST vtx "<<pumi_ment_getGlobalID(e)
             <<" ("<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<")\n";
    else
      std::cout<<"("<<PCU_Comm_Self()<<") vtx "<<pumi_ment_getGlobalID(e)
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
        std::cout<<"("<<PCU_Comm_Self()<<") GHOST e [d "<<d<<", id "<<global_id<<"] down: ";
      else
        std::cout<<"("<<PCU_Comm_Self()<<") e [d "<<d<<", id "<<global_id<<"] down: ";
      for (int i=0; i<num_down; ++i)
        std::cout<<pumi_ment_getGlobalID(down[i])<<" ";
      std::cout<<"\n";
      print_copies(m,e);
    }
    m->end(eit);
  }
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

  destroy_global_numbering(m);
  m->destroyNative();
  apf::destroyMesh(m);
}

void pumi_mesh_verify(pMesh m)
{
  if (pumi::instance()->ghosted_tag)
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" not supported with ghosted mesh\n";
    return;
  }
  apf::verify(m);
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
  assert(parts_vec[i].size()>0);
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
      Up upward;
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
