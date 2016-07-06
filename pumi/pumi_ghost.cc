/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <iostream>
#include <vector>
#include <mpi.h>
#include <PCU.h>
#include <map>
#include <set>
#include <assert.h>
#include <malloc.h>

#include "apf.h"
#include "apfMDS.h"

using std::map;
using std::set;

Ghosting::Ghosting(pMesh mesh, int d)
{
  m = mesh;
  ghost_dim = d;
  ghosted_tag = m->createIntTag("_ghosted_",1);
  ghost_tag = m->createIntTag("_ghost_",1);
}

Ghosting::~Ghosting()
{
/*  map<pMeshEnt, set<int> >::iterator mapit;
  for (mapit=pid_map->begin(); mapit!= pid_map->end(); ++mapit)
    mesh->removeTag(mapit->first,tag);
*/
  for (int i=0; i<4; ++i)
    pid_map[i].clear();
}

int Ghosting::count(pMeshEnt e)
{
  int ent_dim = apf::getDimension(m, e);
  if (has(e))
    pid_map[ent_dim].size();
  else
    return 0;
}

int Ghosting::count(int d)
{
  return pid_map[d].size();
}

pMeshEnt Ghosting::get(int d, int i)
{
  assert(i<count(d));
  map<pMeshEnt, Parts >::iterator mapit = pid_map[d].begin();
  for (int k=0; k<i; ++k)
    mapit++;
  return mapit->first;
}

bool Ghosting::has(pMeshEnt e)
{
  int ent_dim = apf::getDimension(m, e);
  map<pMeshEnt, Parts >::iterator mapit = pid_map[ent_dim].find(e);
  if (mapit!=pid_map[ent_dim].end())
    return true;
  else
    return false;
}

void Ghosting::send(pMeshEnt e, int to)
{
  int ent_dim = apf::getDimension(m, e);
  assert(ent_dim==ghost_dim && to!=PCU_Comm_Self());
  if (!m->hasTag(e,ghost_tag)) 
    pid_map[ent_dim][e].insert(to);
}

/** assign a destination part id of all entities of dimension */
void Ghosting::send(int d, int to)
{
  pMeshEnt e;
  assert(to!=PCU_Comm_Self() && d==ghost_dim);
  apf::MeshIterator* it = m->begin(d);
  while ((e = m->iterate(it)))
    pid_map[d][e].insert(to);
  m->end(it);
}

// *****************************************
void copy_pids(Ghosting* plan, pMeshEnt s, pMeshEnt d)
// *****************************************
{ 
  int s_dim = apf::getDimension(plan->getMesh(), s);
  int d_dim = apf::getDimension(plan->getMesh(), d);
  APF_ITERATE(Parts, plan->pid_map[s_dim][s], pit)
    plan->pid_map[d_dim][d].insert(*pit);
}

static void getAffected(
    pMesh m,
    Ghosting* plan,
    EntityVector affected[4])
{
  int maxDimension = plan->ghost_dim;
  int self = PCU_Comm_Self();
  affected[maxDimension].reserve(plan->count(maxDimension));
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
  for (mapit=plan->pid_map[plan->ghost_dim].begin();mapit!=plan->pid_map[plan->ghost_dim].end();++mapit)
    affected[maxDimension].push_back(mapit->first);

  int dummy=1;
  pTag tag = m->createIntTag("ghosting_affected",1);
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
          // copy destination pids
          copy_pids(plan, up, adjacent[i]);
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

// packing and unpacking are defined in apfMigrate.cc

// *****************************************
void unify_pids(Ghosting* plan, EntityVector affected[4])
// *****************************************
{
  pMeshEnt e;
  int self = PCU_Comm_Self();

  void* msg_send;
  pMeshEnt* s_ent;
  size_t msg_size;

  for (int d=0; d <plan->ghost_dim;++d)
  {
    PCU_Comm_Begin();
    APF_ITERATE(EntityVector,affected[d],it)
    {
      e = *it;
      apf::Copies remotes;
      plan->getMesh()->getRemotes(e,remotes);
      int num_pids=plan->pid_map[d][e].size();
      APF_ITERATE(apf::Copies,remotes,rit)
      {
        msg_size=sizeof(pMeshEnt) +num_pids*sizeof(int);
        msg_send = malloc(msg_size);
        
        s_ent = (pMeshEnt*)msg_send; 
        *s_ent = rit->second; 

        int *pids = (int*)((char*)msg_send + sizeof(pMeshEnt));
        int i = 0;
        APF_ITERATE(Parts, plan->pid_map[d][e], pit)
        {
          pids[i]=*pit;
          ++i;
        }
        PCU_Comm_Write(rit->first, (void*)msg_send, msg_size);
        free(msg_send); 
      }
    }  
    PCU_Comm_Send();
    // receive phase
    void *msg_recv;
    int pid_from, r_dim;
    int* pids;
    pMeshEnt r;
    while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
    {
      r = *((pMeshEnt*)msg_recv); 
      r_dim = apf::getDimension(plan->getMesh(), r);
      pids = (int*)((char*)msg_recv+sizeof(pMeshEnt)); 
      int num_pids = (msg_size-sizeof(pMeshEnt))/sizeof(int);
      for (int i = 0; i < num_pids; ++i)
      {
        if (pids[i]==pumi_rank()) continue;
        plan->pid_map[r_dim][r].insert(pids[i]);
        std::cout<<"("<<pumi_rank()<<") "<<__func__<<": plan->pid_map["<<r_dim<<"]["<<pumi_ment_getglobalid(r)<<"].insert("<<pids[i]<<");\n";
      }
    }
  }//dimension loop
}

static pMeshEnt unpackGhost(Ghosting* plan, apf::DynamicArray<pTag>& tags)
{
  int from = PCU_Comm_Sender();
  int type;
  PCU_COMM_UNPACK(type);
  pMeshEnt sender;
  apf::ModelEntity* c;
  Parts residence;
  apf::unpackCommon(plan->getMesh(),sender,c,residence);
  pMeshEnt entity;
  if (type == apf::Mesh::VERTEX)
    entity = apf::unpackVertex(plan->getMesh(),c);
  else
    entity = apf::unpackNonVertex(plan->getMesh(),type,c);
  plan->getMesh()->setResidence(entity,residence);
  apf::unpackTags(plan->getMesh(),entity,tags);

  // tag ghoating flag
  int dummy=1;
  plan->getMesh()->setIntTag(entity,plan->ghost_tag,&dummy);
  /* store the sender as a ghost copy */
  plan->getMesh()->addGhost(entity, from, sender);
  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": ghost (from "<<from<<", dim "<<type<<", id "<< pumi_ment_getglobalid(entity)<<") \n";
  return entity;
}

// *********************************************************
static void sendGhosts(Ghosting* plan, int ent_dim,
    EntityVector& senders, apf::DynamicArray<pTag>& tags)
// *********************************************************
{
  int dummy=1;
  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<__LINE__<<", senders.size() "<<senders.size()<<"\n";
  APF_ITERATE(EntityVector,senders,it)
  {
    pMeshEnt e = *it;
    Copies remotes;
    plan->getMesh()->getRemotes(e,remotes);
    Parts ghosting_pids;
    APF_ITERATE(Parts, plan->pid_map[ent_dim][e], pit)
      ghosting_pids.insert(*pit);

    Parts sendTo;
    apf::split(remotes,ghosting_pids,sendTo);
    if (!sendTo.size())
    {
      std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<__LINE__<<", sendTo.size()=0 for ent "<<pumi_ment_getglobalid(e)<<"\n";
      continue;
    }
    APF_ITERATE(Parts,sendTo,sit)
    {
      std::cout<<"("<<pumi_rank()<<") "<<__func__<<": pack dim "<<ent_dim<<", id "<<pumi_ment_getglobalid(e)<<" to "<<*sit<<"\n";
      apf::packEntity(plan->getMesh(),*sit,e,tags);
    }
    // attach ghosted_tag after packing to avoid this tag migrated
    plan->getMesh()->setIntTag(e,plan->ghosted_tag,&dummy);
  }
}

// *********************************************************
static void receiveGhosts(Ghosting* plan, apf::DynamicArray<pTag>& tags,
    EntityVector& received)
// *********************************************************
{
  received.reserve(1024);
  while (PCU_Comm_Receive())
    received.push_back(unpackGhost(plan,tags));
}

// *********************************************************
static void echoGhostCopy(pMesh m, EntityVector& received)
// *********************************************************
{
  APF_ITERATE(EntityVector,received,it)
  {
    pMeshEnt entity = *it;
    /* the remote copies are currently temporary
       storage for the sender */
    apf::Copies temp;
    m->getGhosts(entity,temp);
    int from = temp.begin()->first;
    pMeshEnt sender = temp.begin()->second;
    PCU_COMM_PACK(from,sender);
    PCU_COMM_PACK(from,entity);
    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": echo dim "<<apf::getDimension(m, entity)<<" id "<<pumi_ment_getglobalid(entity)<<" "<<getMdsIndex(m, entity)<<" to "<<from<<"\n";
  }
}

// *********************************************************
static void receiveGhostCopy(pMesh m)
// *********************************************************
{
  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while ( ! PCU_Comm_Unpacked())
    {
      pMeshEnt sender;
      PCU_COMM_UNPACK(sender);
      pMeshEnt entity;
      PCU_COMM_UNPACK(entity);
      assert(entity);
      std::cout<<"("<<pumi_rank()<<") "<<__func__<<": dim "<<apf::getDimension(m, entity)<<" id "<<pumi_ment_getglobalid(entity)<<" "<<getMdsIndex(m, entity)<<" from "<<from<<"\n";
      m->addGhost(sender, from, entity);
    }
  }
}

// *********************************************************
static void setupGhosts(pMesh m, EntityVector& received, EntityVector& senders)
// *********************************************************
{
  PCU_Comm_Begin();
  echoGhostCopy(m,received);
  PCU_Comm_Send();
  receiveGhostCopy(m);
}

// *********************************************************
static void moveGhosts(Ghosting* plan, EntityVector senders[4])
// *********************************************************
{
  apf::DynamicArray<pTag> tags;
  plan->getMesh()->getTags(tags);

  for (int dimension = 0; dimension <= plan->ghost_dim; ++dimension)
  {
    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": START dim "<<dimension<<"\n";
    PCU_Comm_Begin();
    sendGhosts(plan, dimension, senders[dimension], tags);
    PCU_Comm_Send();
    EntityVector received;
    receiveGhosts(plan,tags,received);
    setupGhosts(plan->getMesh(),received,senders[dimension]);
    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": END dim "<<dimension<<"\n";
  }
}


// *********************************************************
void pumi_ghost_create(pMesh m, Ghosting* plan)
// *********************************************************
{
  EntityVector affected[4];
  getAffected(m, plan,affected);
  EntityVector senders[4];
  apf::getSenders(m,affected,senders);
  unify_pids(plan, affected);
  moveGhosts(plan,senders);
//  updateMatching(m,affected,senders);
//  deleteOldEntities(m,affected);
  m->acceptChanges();
  delete plan;
}


// *********************************************************
pMesh pumi_ghost_createlayer (int brgType, int ghostType, int numLayer, int includeCopy)
// *********************************************************
{
  
  if (brgType!=0 && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: bridge type "<<brgType<<" not supported\n";
    return NULL;
  }
  
  if (ghostType!=pumi::instance()->mesh->getDimension() && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" faild: ghost type "<<brgType<<" not supported\n";
    return NULL;
  }
  
  if (numLayer<1 && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid numLayer "<<numLayer<<"\n";
    return NULL;
  }

  if (includeCopy==0 && !pumi_rank()) 
  {
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: includeCopy=0"<<" not supported\n";
    return NULL;
  }
  if (!pumi_rank()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: accumulative ghosting not supported\n";
      return NULL;
}

// *********************************************************
void pumi_ghost_delete (pMesh m)
// *********************************************************
{
  // remove ghost copies
  // remove "ghosts" and "ghosted_tag" from ghosted copy
  pTag ghosted_tag=m->findTag("_ghosted_");
  pTag ghost_tag=m->findTag("_ghost_");
  if (ghost_tag) m->destroyTag(ghost_tag);
  if (ghosted_tag) m->destroyTag(ghosted_tag);

  if (!pumi_rank()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: not supported\n";
  return;
}

// *********************************************************
void pumi_ghost_info (pMesh m, std::vector<int>& ghostinfo)
// *********************************************************
{
  if (!pumi_rank()) 
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: not supported\n";
}

