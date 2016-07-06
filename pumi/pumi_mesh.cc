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
#include <map>
#include <assert.h>
#include <malloc.h>
#include "apf.h"

using std::map;

void generate_globalid(pMesh m, pTag tag, int dim)
{
  pMeshEnt e;
  int own_partid, num_own=0, myrank=PCU_Comm_Self();

  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    own_partid = m->getOwner(e);
    if (own_partid==myrank)
      ++num_own;
  }
  m->end(it);

  PCU_Exscan_Ints(&num_own,1);
  int initial_id=num_own;

  PCU_Comm_Begin();
  it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    own_partid = m->getOwner(e);
    if (own_partid==myrank)
    {
      m->setIntTag(e, tag, &initial_id);
      Copies remotes;
      m->getRemotes(e, remotes);
      APF_ITERATE(Copies, remotes, it)
      {
        PCU_COMM_PACK(it->first, it->second);
        PCU_Comm_Pack(it->first, &initial_id, sizeof(int));
      }
      ++initial_id;
    }
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
void generate_global_numbering(apf::Mesh2* m)
//*******************************************************

{
  pTag tag = m->findTag("global_id");
  if (tag)  // destroy existing tag
  {
    for (int i=0; i<4; ++i)
      apf::removeTagFromDimension(m, tag, m->getDimension());
  }  
  else
    tag = m->createIntTag("global_id",1);

  for (int i=0; i<4; ++i)
    generate_globalid(m, tag, i);
}

//*******************************************************
void destroy_global_numbering(apf::Mesh2* m)
//*******************************************************
{
  pTag tag = m->findTag("global_id");
  for (int i=0; i<4; ++i)
    apf::removeTagFromDimension(m, tag, m->getDimension());

  m->destroyTag(tag);
}

// *********************************************************
pumi::pumi(): mesh(NULL), model(NULL) {}

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
  pTag weights = Parma_WeighByMemory(m);
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



// load a serial mesh on master process then distribute as per the distribution object
pMesh pumi_mesh_loadserial(pGeom g, const char* filename, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  // set proc_group_id
  int self = PCU_Comm_Self();

  MPI_Comm prevComm = PCU_Get_Comm();
  int num_target_part = PCU_Comm_Peers();
  bool isMaster = ((PCU_Comm_Self() % num_target_part) == 0);
  pMesh m = 0;
  split_comm(num_target_part);
  if (isMaster) 
    m = apf::loadMdsMesh(g, filename);
  merge_comm(prevComm);
  pumi::instance()->mesh = expandMdsMesh(m, g, 1);
  generate_global_numbering(pumi::instance()->mesh);
  return pumi::instance()->mesh;
}


pMesh pumi_mesh_load(pGeom g, const char* filename, int num_in_part, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  // set proc_group_id
  int self = PCU_Comm_Self();

  if (num_in_part==1) // do static partitioning
  {
    MPI_Comm prevComm = PCU_Get_Comm();
    int num_target_part = PCU_Comm_Peers()/num_in_part;
    bool isMaster = ((PCU_Comm_Self() % num_target_part) == 0);
    pMesh m = 0;
    apf::Migration* plan = 0;   
    split_comm(num_target_part);
    if (isMaster) {
      m = apf::loadMdsMesh(g, filename);
      plan = getPlan(m, num_target_part);
    }
    merge_comm(prevComm);
    pumi::instance()->mesh = apf::repeatMdsMesh(m, g, plan, num_target_part);
  }
  else
    pumi::instance()->mesh = apf::loadMdsMesh(g, filename);

  generate_global_numbering(pumi::instance()->mesh);
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
  if (!PCU_Comm_Self()) std::cout<<"\n=== mesh size and tag info === \nglobal ";
  printStats(m);

  int* local_entity_count = new int[4*PCU_Comm_Peers()];
  for (int i=0; i<4*PCU_Comm_Peers();++i)
    local_entity_count[i]=0;
  for (int d=0; d<4;++d)
    local_entity_count[4*pumi_rank()+d] = m->count(d);
  
  int* global_entity_count = new int[4*PCU_Comm_Peers()]; 

  MPI_Allreduce(local_entity_count, global_entity_count, 4*PCU_Comm_Peers(), MPI_INT, MPI_SUM, PCU_Get_Comm());
  if (pumi_rank()) return;

  for (int p=0; p<PCU_Comm_Peers(); ++p)
    std::cout<<"(p"<<p<<") # local ent: v "<<global_entity_count[p*4]
        <<", e "<<global_entity_count[p*4+1]
        <<", f "<<global_entity_count[p*4+2]
        <<", r "<<global_entity_count[p*4+3]<<"\n";

  delete [] local_entity_count;
  delete [] global_entity_count;

  apf::DynamicArray<pTag> tags;
  m->getTags(tags);
  int n = tags.getSize();
  for (int i = 0; i < n; ++i) 
    std::cout<<"tag "<<i<<": \""<< m->getTagName(tags[i])<<"\", type "<< m->getTagType(tags[i])<<", size "<< m->getTagSize(tags[i])<<"\n";
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
  destroy_global_numbering(m);
  m->destroyNative();
  apf::destroyMesh(m);
}

void pumi_mesh_verify(pMesh m)
{
  apf::verify(m);
}

Distribution::Distribution(pMesh m)
{
  mesh = m;
}

int Distribution::count()
{
  return element_map.size();
}

pMeshEnt Distribution::get(int i)
{
  map<pMeshEnt, Parts >::iterator mapit = element_map.begin();
  for (int k=0; k<i; ++k)
    mapit++;
  return mapit->first;
}

bool Distribution::has(pMeshEnt e)
{
  map<pMeshEnt, Parts >::iterator mapit = element_map.find(e);
  if (mapit!=element_map.end())
    return true;
  else
    return false;
}

void Distribution::send(pMeshEnt e, int to)
{
  if (to!=PCU_Comm_Self())
    element_map[e].insert(to);
}

Parts& Distribution::sending(pMeshEnt e)
{
  return element_map[e];
}

void Distribution::print()
{
  map<pMeshEnt, Parts >::iterator mapit;
  for (mapit = element_map.begin(); mapit!=element_map.end(); ++mapit)
  {
    int ent_id = pumi_ment_getglobalid(mapit->first);
    APF_ITERATE(Parts,element_map[mapit->first],pit)
      std::cout<<"("<<PCU_Comm_Self()<<") send e "<<ent_id<<" to "<<*pit<<"\n";
  }
}

static void distr_getAffected (pMesh m, Distribution* plan, EntityVector affected[4])
{
  int maxDimension = m->getDimension();
  int self = PCU_Comm_Self();
  affected[maxDimension].reserve(plan->count());
/*
  for (int i=0; i < plan->count(); ++i)
  {
    pMeshEnt e = plan->get(i);
    if (plan->sending(e) != self) {
      assert(apf::getDimension(m, e) == m->getDimension());
      affected[maxDimension].push_back(e);
    }
  }
*/
  map<pMeshEnt, Parts >::iterator mapit;
  for (mapit=plan->element_map.begin();mapit!=plan->element_map.end();++mapit)
    affected[maxDimension].push_back(mapit->first);

  int dummy=1;
  pTag tag = m->createIntTag("distribution_affected",1);
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
  for (int i=0; i < plan->count(); ++i)
  {
    pMeshEnt e = plan->get(i);
    Parts res = distr_makeResidence(plan->sending(e));
    m->setResidence(e,res);
  }
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
pMesh distr_repeatMdsMesh(pMesh m, pGeom g, Distribution* plan,  int factor)
// *********************************************************
{
  double t0 = PCU_Time();
  if (PCU_Comm_Self() % factor != 0)
    plan = new Distribution(m);
  distribute(m, plan);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("[PUMI INFO] mesh distributed from %d to %d in %f seconds\n",
        PCU_Comm_Peers() / factor, PCU_Comm_Peers(), t1 - t0);
  return m;
}

// *********************************************************
void pumi_mesh_distribute(pMesh m, Distribution* plan)
// *********************************************************
{
  if (PCU_Comm_Peers()==1) return;
  distribute(m, plan);
}
