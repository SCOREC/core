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
  map<pMeshEnt, set<int> >::iterator mapit = pid_map[d].begin();
  for (int k=0; k<i; ++k)
    mapit++;
  return mapit->first;
}

bool Ghosting::has(pMeshEnt e)
{
  int ent_dim = apf::getDimension(m, e);
  map<pMeshEnt, set<int> >::iterator mapit = pid_map[ent_dim].find(e);
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
  std::set<int>::iterator pidit;
  for (pidit=plan->pid_map[s_dim][s].begin(); pidit!=plan->pid_map[s_dim][s].end(); ++pidit)
    plan->pid_map[d_dim][d].insert(*pidit);
}

typedef std::vector<pMeshEnt> EntityVec;

// *****************************************
static void getAffected(Ghosting* plan, EntityVec affected[4])
// *****************************************
{  
  int self = PCU_Comm_Self();
  affected[plan->ghost_dim].reserve(plan->count(plan->ghost_dim));
  map<pMeshEnt, set<int> >::iterator mapit;
  for (mapit=plan->pid_map[plan->ghost_dim].begin();mapit!=plan->pid_map[plan->ghost_dim].end();++mapit)
    affected[plan->ghost_dim].push_back(mapit->first);

  int dummy=1;
  pTag tag = plan->getMesh()->createIntTag("apf_migrate_affected",1);
  for (int dimension=plan->ghost_dim-1; dimension >= 0; --dimension)
  {
    int upDimension = dimension + 1;
    PCU_Comm_Begin();
    APF_ITERATE(EntityVec,affected[upDimension],it)
    {
      pMeshEnt up = *it;
      apf::Downward adjacent;
      int na = plan->getMesh()->getDownward(up,dimension,adjacent);
      for (int i=0; i < na; ++i)
      {
        if (!plan->getMesh()->hasTag(adjacent[i],tag))
        {
          plan->getMesh()->setIntTag(adjacent[i],tag,&dummy);
          affected[dimension].push_back(adjacent[i]);
          copy_pids(plan, up, adjacent[i]);
        }
        apf::Copies remotes;
        plan->getMesh()->getRemotes(adjacent[i],remotes);
        APF_ITERATE(apf::Copies,remotes,rit)
          PCU_COMM_PACK(rit->first,rit->second);
        if (plan->getMesh()->hasMatching())
        {
          apf::Matches matches;
          plan->getMesh()->getMatches(adjacent[i],matches);
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
      if (!plan->getMesh()->hasTag(entity,tag))
      {
        plan->getMesh()->setIntTag(entity,tag,&dummy);
        affected[dimension].push_back(entity);
      }
    }
    APF_ITERATE(EntityVec,affected[dimension],it)
      plan->getMesh()->removeTag(*it,tag);
  }//dimension loop
  plan->getMesh()->destroyTag(tag);
}


static void getSenders(
    pMesh m,
    EntityVec affected[4],
    EntityVec senders[4])
{
  for (int i=0; i < 4; ++i)
  {
    /* maybe overkill for pre-allocation, but
       ensures the vector will never re-allocate. */
    senders[i].reserve(affected[i].size());
    APF_ITERATE(EntityVec,affected[i],it)
      if (m->isOwned(*it))
        senders[i].push_back(*it);
  }
}

/* at this point if one matched copy is affected, all of
   its matches are in the affected set and the "senders", or
   representatives of a set of remote copies, are in
   the senders vector.
   So, if each sender just remembers the other senders
   it is matched to, that is enough for them to re-negotiate
   matches after each of them does the job of creating new
   entities and computing new remote copies */
static void reduceMatchingToSenders(
    pMesh m,
    EntityVec senders[4])
{
  if ( ! m->hasMatching()) return;
  PCU_Comm_Begin();
  for (int d=0; d < 4; ++d)
  {
    for (size_t i=0; i < senders[d].size(); ++i)
    {
      pMeshEnt e = senders[d][i];
      apf::Matches matches;
      m->getMatches(e,matches);
      for (size_t j=0; j < matches.getSize(); ++j)
      { /* advertise to the match that this is a sender */
        int to = matches[j].peer;
        PCU_COMM_PACK(to,matches[j].entity);
        PCU_COMM_PACK(to,e);
      }
      /* now forget all other matchings */
      m->clearMatches(e);
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    int sender = PCU_Comm_Sender();
    while ( ! PCU_Comm_Unpacked())
    {
      pMeshEnt e;
      PCU_COMM_UNPACK(e);
      pMeshEnt match;
      PCU_COMM_UNPACK(match);
/* all matches of the sender receive this.
   We only care that the senders themselves receive
   messages from only other senders matched to them,
   but some non-senders do harmlessly add a few matchings */
      m->addMatch(e,sender,match);
    }
  }
/* now senders only match other senders */
}

// packing

static void packParts(int to, apf::Parts& parts)
{
  size_t n = parts.size();
  PCU_COMM_PACK(to,n);
  APF_ITERATE(apf::Parts,parts,it)
  {
    int p = *it;
    PCU_COMM_PACK(to,p);
  }
}

static void packCommon(
    pMesh m,
    int to,
    pMeshEnt e)
{
  PCU_COMM_PACK(to,e);
  apf::ModelEntity* me = m->toModel(e);
  int modelType = m->getModelType(me);
  PCU_COMM_PACK(to,modelType);
  int modelTag = m->getModelTag(me);
  PCU_COMM_PACK(to,modelTag);
  apf::Parts residence;
  m->getResidence(e,residence);
  packParts(to,residence);
}

static void packVertex(
    pMesh m,
    int to,
    pMeshEnt e)
{
  apf::Vector3 p;
  m->getPoint(e,0,p);
  PCU_COMM_PACK(to,p);
  m->getParam(e,p);
  PCU_COMM_PACK(to,p);
}

static void packReference(
    pMesh m,
    int to,
    pMeshEnt e)
{
  apf::Copies remotes;
  m->getRemotes(e,remotes);
  apf::Copies::iterator found = remotes.find(to);
  pMeshEnt remote = found->second;
  assert(remote);
  PCU_COMM_PACK(to,remote);
}

static void packDownward(pMesh m, int to, pMeshEnt e)
{
  apf::Downward down;
  int d = getDimension(m, e);
  int n = m->getDownward(e,d-1,down);
  PCU_COMM_PACK(to,n);
  for (int i=0; i < n; ++i)
    packReference(m,to,down[i]);
}

static void packNonVertex(
    pMesh m,
    int to,
    pMeshEnt e)
{
  packDownward(m,to,e);
}

static void packTags(
    pMesh m,
    int to,
    pMeshEnt e,
    apf::DynamicArray<pTag>& tags)
{
  size_t total = tags.getSize();
  size_t n = 0;
  for (size_t i=0; i < total; ++i)
    if (m->hasTag(e,tags[i]))
      ++n;
  PCU_COMM_PACK(to,n);
  for (size_t i=0; i < total; ++i)
  {
    pTag tag = tags[i];
    if (m->hasTag(e,tag))
    {
      PCU_COMM_PACK(to,i);
      int type = m->getTagType(tag);
      int size = m->getTagSize(tag);
      if (type == apf::Mesh2::DOUBLE)
      {
        apf::DynamicArray<double> d(size);
        m->getDoubleTag(e,tag,&(d[0]));
        PCU_Comm_Pack(to,&(d[0]),size*sizeof(double));
      }
      if (type == apf::Mesh2::INT)
      {
        apf::DynamicArray<int> d(size);
        m->getIntTag(e,tag,&(d[0]));
        PCU_Comm_Pack(to,&(d[0]),size*sizeof(int));
      }
    }
  }
}

// unpacking
static void unpackParts(apf::Parts& parts)
{
  size_t n;
  PCU_COMM_UNPACK(n);
  for (size_t i=0;i<n;++i)
  {
    int p;
    PCU_COMM_UNPACK(p);
    parts.insert(p);
  }
}

static void unpackCommon(
    pMesh m,
    pMeshEnt& sender,
    apf::ModelEntity*& c,
    apf::Parts& residence)
{
  PCU_COMM_UNPACK(sender);
  int modelType,modelTag;
  PCU_COMM_UNPACK(modelType);
  PCU_COMM_UNPACK(modelTag);
  c = m->findModelEntity(modelType,modelTag);
  unpackParts(residence);
}

static pMeshEnt unpackVertex(
    pMesh m,
    apf::ModelEntity* c)
{
  apf::Vector3 point;
  PCU_COMM_UNPACK(point);
  apf::Vector3 param;
  PCU_COMM_UNPACK(param);
  return m->createVertex(c,point,param);
}

static void unpackDownward(
    apf::Downward& entities)
{
  int n;
  PCU_COMM_UNPACK(n);
  for (int i=0; i < n; ++i)
    PCU_COMM_UNPACK(entities[i]);
}

static pMeshEnt unpackNonVertex(
    pMesh m,
    int type, apf::ModelEntity* c)
{
  apf::Downward down;
  unpackDownward(down);
  return m->createEntity(type,c,down);
}


static void unpackTags(
    pMesh m,
    pMeshEnt e,
    apf::DynamicArray<pTag>& tags)
{
  size_t n;
  PCU_COMM_UNPACK(n);
  for (size_t t=0; t < n; ++t)
  {
    size_t i;
    PCU_COMM_UNPACK(i);
    pTag tag = tags[i];
    int type = m->getTagType(tag);
    int size = m->getTagSize(tag);
    if (type == apf::Mesh2::DOUBLE)
    {
      apf::DynamicArray<double> d(size);
      PCU_Comm_Unpack(&(d[0]),size*sizeof(double));
      m->setDoubleTag(e,tag,&(d[0]));
    }
    if (type == apf::Mesh2::INT)
    {
      apf::DynamicArray<int> d(size);
      PCU_Comm_Unpack(&(d[0]),size*sizeof(int));
      m->setIntTag(e,tag,&(d[0]));
    }
  }
}


// *****************************************
void update_pids(Ghosting* plan, EntityVec affected[4])
// *****************************************
{
  pMeshEnt e;
  int self = PCU_Comm_Self();
  set<int>::iterator setit;

  void* msg_send;
  pMeshEnt* s_ent;
  size_t msg_size;

  for (int d=0; d <plan->ghost_dim;++d)
  {
    PCU_Comm_Begin();
    APF_ITERATE(EntityVec,affected[d],it)
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
        for (setit=plan->pid_map[d][e].begin();setit!=plan->pid_map[d][e].end();++setit)
        {
          pids[i]=*setit;
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
        plan->pid_map[r_dim][r].insert(pids[i]);
    }
  }//dimension loop
}

static void packGhost(
    pMesh m,
    int to,
    pMeshEnt e,
    apf::DynamicArray<pTag>& tags)
{
  int type = m->getType(e);
  PCU_COMM_PACK(to,type);
  packCommon(m,to,e);
  if (type == apf::Mesh::VERTEX)
    packVertex(m,to,e);
  else
    packNonVertex(m,to,e);
  packTags(m,to,e,tags);
}

static pMeshEnt unpackGhost(
    Ghosting* plan,
    apf::DynamicArray<pTag>& tags)
{
  int from = PCU_Comm_Sender();
  int type, dummy=1;
  PCU_COMM_UNPACK(type);
  pMeshEnt sender;
  apf::ModelEntity* c;
  apf::Parts residence;
  unpackCommon(plan->getMesh(),sender,c,residence);
  pMeshEnt entity;
  if (type == apf::Mesh::VERTEX)
    entity = unpackVertex(plan->getMesh(),c);
  else
    entity = unpackNonVertex(plan->getMesh(),type,c);
  plan->getMesh()->setIntTag(entity,plan->ghost_tag,&dummy);
  plan->getMesh()->setResidence(entity,residence);
  unpackTags(plan->getMesh(),entity,tags);

  /* store the sender as a ghost copy */
  plan->getMesh()->addGhost(entity, from, sender);
  return entity;
}


static void sendGhosts(
    Ghosting* plan,
    int ent_dim,
    EntityVec& senders,
    apf::DynamicArray<pTag>& tags)
{
  int dummy=1;
  APF_ITERATE(EntityVec,senders,it)
  {
    pMeshEnt entity = *it;
    apf::Copies remotes;
    plan->getMesh()->getRemotes(entity,remotes);
    apf::Parts residence;
    plan->getMesh()->getResidence(entity,residence);
    apf::Parts sendTo;
    set<int>::iterator setit;
    APF_ITERATE(apf::Parts,residence,pit)
    {
      setit=plan->pid_map[ent_dim][entity].find(*pit);
      if (setit==plan->pid_map[ent_dim][entity].end())
        sendTo.insert(*pit);
    }
    APF_ITERATE(apf::Parts,sendTo,sit)
      packGhost(plan->getMesh(),*sit,entity,tags);
    // attach ghosted_tag after packing to avoid this tag migrated
    plan->getMesh()->setIntTag(entity,plan->ghosted_tag,&dummy);
  }
}

// *********************************************************
static void receiveGhosts(Ghosting* plan, apf::DynamicArray<pTag>& tags,
    EntityVec& received)
// *********************************************************
{
  received.reserve(1024);
  while (PCU_Comm_Receive())
    received.push_back(unpackGhost(plan,tags));
}

// *********************************************************
static void _echo_(pMesh m, EntityVec& received)
// *********************************************************
{
  APF_ITERATE(EntityVec,received,it)
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
  }
}

// *********************************************************
static void _receive_(pMesh m)
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
      m->addGhost(sender, from, entity);
    }
  }
}

// *********************************************************
static void setupGhosts(pMesh m, EntityVec& received, EntityVec& senders)
// *********************************************************
{
  PCU_Comm_Begin();
  _echo_(m,received);
  PCU_Comm_Send();
  _receive_(m);
}

// *********************************************************
static void moveGhosts(Ghosting* plan, EntityVec senders[4])
// *********************************************************
{
  apf::DynamicArray<pTag> tags;
  plan->getMesh()->getTags(tags);

  for (int dimension = 0; dimension <= plan->ghost_dim; ++dimension)
  {
    PCU_Comm_Begin();
    sendGhosts(plan, dimension, senders[dimension], tags);
    PCU_Comm_Send();
    EntityVec received;
    receiveGhosts(plan,tags,received);
    setupGhosts(plan->getMesh(),received,senders[dimension]);
  }
}


// *********************************************************
void pumi_ghost_create(pMesh m, Ghosting* plan)
// *********************************************************
{
  EntityVec affected[4];
  getAffected(plan,affected);
  EntityVec senders[4];
  getSenders(m,affected,senders);
//  reduceMatchingToSenders(m,senders);
//  updateResidences(m,plan,affected);
  update_pids(plan, affected);
  delete plan;
  moveGhosts(plan,senders);
//  updateMatching(m,affected,senders);
//  deleteOldEntities(m,affected);
  m->acceptChanges();
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

