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

  parts_index_tag = m->findTag("_parts_index_");
  if (!parts_index_tag) 
    parts_index_tag = m->createIntTag("_parts_index_", 1);
  else // clean the existing parts_index_tag;
  {
    for (int d=0; d<4; ++d)
      apf::removeTagFromDimension(m, parts_index_tag, d);
  }
}

Ghosting::~Ghosting()
{
  for (int i=0; i<4; ++i)
  {
    for (std::vector<Parts*>::iterator vit=parts_vec[i].begin(); vit!=parts_vec[i].begin(); ++vit)
      delete *vit;
    parts_vec[i].clear();
  }
  // FIXME: delete tag
  parts_index_tag = m->findTag("_parts_index_");
  if (parts_index_tag) 
  {
    if (!PCU_Comm_Self()) std::cout<<"deleting parts_vec_index from entities\n";
    for (int d=0; d<=ghost_dim; ++d)
      apf::removeTagFromDimension(m, parts_index_tag, d);
  }
}

bool Ghosting::has(pMeshEnt e)
{
  if (m->hasTag(e, parts_index_tag))
    return true;
  else
    return false;
}

void Ghosting::send(pMeshEnt e, int to)
{
  if (to==PCU_Comm_Self()) return;

  int d = getDimension(m, e);
  int index=-1;
  if (!m->hasTag(e, parts_index_tag))
  {
    index=parts_vec[d].size();
    m->setIntTag(e,parts_index_tag,&index);    
    parts_vec[d].push_back(new Parts);
  }
  else
    m->getIntTag(e, parts_index_tag,&index);
  assert(index!=-1);
//  std::cout<<"("<<PCU_Comm_Self()<<") send e (dim "<<d<<", id "<<pumi_ment_getglobalid(e)<<") to "<<to<<")\n"; 
  parts_vec[d][index]->insert(to);
}

/** assign a destination part id of all entities of dimension */
void Ghosting::send (int to)
{
  if (to==PCU_Comm_Self()) return;

  pMeshEnt e;
  apf::MeshIterator* it = m->begin(ghost_dim);
  while ((e = m->iterate(it)))
    send(e, to);
  m->end(it);
}

void Ghosting::print()
{
  pMeshEnt e;
  apf::MeshIterator* it = m->begin(ghost_dim);
  int index;
  while ((e = m->iterate(it)))
  {
    if (!m->hasTag(e, parts_index_tag)) continue;
    m->getIntTag(e, parts_index_tag, &index);
    
    APF_ITERATE(Parts,*(parts_vec[ghost_dim][index]),pit)
      std::cout<<"("<<PCU_Comm_Self()<<") ghost e "<<pumi_ment_getglobalid(e)<<" to "<<*pit<<"\n";
  }
  m->end(it);
}

Parts& Ghosting::sending(pMeshEnt e, int d)
{
  int index;
  if(!m->hasTag(e, parts_index_tag))
    std::cout<<"("<<PCU_Comm_Self()<<") ERROR: ghost e (dim "<<d<<") "<<pumi_ment_getglobalid(e)<<" has no parts_index_tag\n";
  assert(m->hasTag(e, parts_index_tag));
  m->getIntTag(e, parts_index_tag, &index);
  return *(parts_vec[d][index]);
}

int Ghosting::count(pMeshEnt e, int d)
{
  if (!m->hasTag(e, parts_index_tag)) return 0;
  int index;
  m->getIntTag(e, parts_index_tag, &index);
  return parts_vec[d][index]->size();
}

int Ghosting::count()
{
  return parts_vec[ghost_dim].size();
}

// *****************************************
static void ghost_getAffected (pMesh m, Ghosting* plan, 
                               EntityVector affected[4])
// *****************************************
{
  int maxDimension = plan->ghost_dim;
  int self = PCU_Comm_Self();
  affected[maxDimension].reserve(plan->count());

  pMeshEnt e;
  apf::MeshIterator* it = m->begin(maxDimension);
  while ((e = m->iterate(it)))
  {
    if (plan->has(e))
      affected[maxDimension].push_back(e);
  }
  m->end(it);

  int dummy=1;
  pTag tag = m->createIntTag("ghosting_affected",1);
  for (int dimension=maxDimension-1; dimension >= 0; --dimension)
  {
    int upDimension = dimension + 1;
    PCU_Comm_Begin();
    APF_ITERATE(EntityVector,affected[upDimension],it)
    {
      pMeshEnt up = *it;
//      std::cout<<"("<<pumi_rank()<<") "<<__func__<<":dim="<<dimension<<", upDim="<<upDimension
//               <<", *up="<<pumi_ment_getglobalid(up)<<"\n";
      apf::Downward adjacent;
      int na = m->getDownward(up,dimension,adjacent);
      for (int i=0; i < na; ++i)
      {
        if ( ! m->hasTag(adjacent[i],tag))
        {
          m->setIntTag(adjacent[i],tag,&dummy);
          affected[dimension].push_back(adjacent[i]);
 //          std::cout<<"("<<pumi_rank()<<") "<<__func__<<": affected["<<dimension<<"].push_back("<<pumi_ment_getglobalid(adjacent[i])<<")\n";
        }
        // copy destination pids
        if (plan->has(up))
          APF_ITERATE(Parts, plan->sending(up, upDimension), pit)
            plan->send(adjacent[i], *pit);//parts_vec[dimension][getMdsIndex(m, adjacent[i])].insert(*pit);
        
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
void ghost_unifyPids(Ghosting* plan, EntityVector affected[4])
// *****************************************
{
  pMesh m = plan->getMesh();
  pMeshEnt e;
  int id, self = PCU_Comm_Self();

  void* msg_send;
  pMeshEnt* s_ent;
  size_t msg_size;
  
  int maxDimension = (plan->ghost_dim==m->getDimension())?plan->ghost_dim-1:plan->ghost_dim;
  for (int d=0; d <=maxDimension;++d)
  {
    PCU_Comm_Begin();
    APF_ITERATE(EntityVector,affected[d],it)
    {
      e = *it;
      if (!plan->has(e)) continue;
      apf::Copies remotes;
      m->getRemotes(e,remotes);
      id = getMdsIndex(m, e);
      int num_pids=plan->count(e, d);
      APF_ITERATE(apf::Copies,remotes,rit)
      {
        msg_size=sizeof(pMeshEnt) +num_pids*sizeof(int);
        msg_send = malloc(msg_size);
        
        s_ent = (pMeshEnt*)msg_send; 
        *s_ent = rit->second; 

        int *pids = (int*)((char*)msg_send + sizeof(pMeshEnt));
        int i = 0;
        APF_ITERATE(Parts, plan->sending(e, d), pit)
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
    int pid_from, r_dim, r_id;
    int* pids;
    pMeshEnt r;
    while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
    {
      r = *((pMeshEnt*)msg_recv); 
      pids = (int*)((char*)msg_recv+sizeof(pMeshEnt)); 
      int num_pids = (msg_size-sizeof(pMeshEnt))/sizeof(int);
      for (int i = 0; i < num_pids; ++i)
        plan->send(r, pids[i]); //parts_vec[r_dim][r_id].insert(pids[i]);
    }
  }//dimension loop
}

// *********************************************************
static pMeshEnt unpackGhost(Ghosting* plan, pTag global_id_tag, apf::DynamicArray<pTag>& tags)
// *********************************************************
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
  // FIXME: how to handle residence for ghost copy?
  plan->getMesh()->setResidence(entity,residence);
  apf::unpackTags(plan->getMesh(),entity,tags);

  // tag ghosting flag
  // attach the global id
  int global_id;
  PCU_Comm_Unpack(&global_id, sizeof(int));
  plan->getMesh()->setIntTag(entity, global_id_tag, &global_id);
  /* store the sender as a ghost copy */
  plan->getMesh()->addGhost(entity, from, sender);
  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<entity<<"(d "<<getDimension(plan->getMesh(), entity)
           <<", id "<<global_id<<")->addGhost("<<from<<", "<<sender<<")\n";
  return entity;
}

// *********************************************************
static void sendGhosts(Ghosting* plan, int ent_dim,
    EntityVector& senders, apf::DynamicArray<pTag>& tags)
// *********************************************************
{
  int dummy=1;
  APF_ITERATE(EntityVector,senders,it)
  {
    pMeshEnt e = *it;
    Copies remotes;
    plan->getMesh()->getRemotes(e,remotes);
    Parts ghosting_pids;
    APF_ITERATE(Parts, plan->sending(e, ent_dim), pit)
      ghosting_pids.insert(*pit);

    Parts sendTo;
    apf::split(remotes,ghosting_pids,sendTo);
    if (!sendTo.size())
      continue;
    int global_id=pumi_ment_getglobalid(e);
    APF_ITERATE(Parts,sendTo,sit)
    {
      
      apf::packEntity(plan->getMesh(),*sit,e,tags);
      PCU_Comm_Pack(*sit,&global_id, sizeof(int));
    }
    // attach ghosted_tag after packing to avoid this tag migrated
//    plan->getMesh()->setIntTag(e,plan->ghosted_tag,&dummy);
  }
}

// *********************************************************
static void receiveGhosts(Ghosting* plan, pTag global_id_tag, apf::DynamicArray<pTag>& tags,
    EntityVector& received)
// *********************************************************
{
  received.reserve(1024);
  while (PCU_Comm_Receive())
    received.push_back(unpackGhost(plan,global_id_tag, tags));
}

// *********************************************************
static void setupGhosts(pMesh m, EntityVector& received, EntityVector& senders)
// *********************************************************
{
  PCU_Comm_Begin();
  APF_ITERATE(EntityVector,received,it)
  {
    pMeshEnt entity = *it;
    /* the remote copies are currently temporary
       storage for the sender */
    apf::Copies temp;
    m->getGhosts(entity,temp);
    int to = temp.begin()->first;
    PCU_COMM_PACK(to,temp.begin()->second); // sender
    PCU_COMM_PACK(to,entity);
//    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": echo dim "<<apf::getDimension(m, entity)<<" id "<<getMdsIndex(m, entity)<<" to "<<to<<", entity="<<entity<<", sender="<<temp.begin()->second<<"\n";
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
  {
    int from = PCU_Comm_Sender();
    pMeshEnt entity;
    PCU_COMM_UNPACK(entity);
    pMeshEnt sender;
    PCU_COMM_UNPACK(sender);
//    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": received entity="<<entity<<", sender="<<sender<<"\n";
  //  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": received dim "<<apf::getDimension(m, entity)<<" id "<<pumi_ment_getglobalid(entity)<<" from "<<from<<"\n";
    m->addGhost(entity, from, sender);
    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<entity<<"(d "<<getDimension(m,entity)
             <<", id "<<pumi_ment_getglobalid(entity)<<")->addGhost("<<from<<", "<<sender<<")\n";
  }
}

// *********************************************************
static void ghost_moveEntities(Ghosting* plan, EntityVector senders[4])
// *********************************************************
{
  apf::DynamicArray<pTag> tags;
  plan->getMesh()->getTags(tags);
  pTag global_id_tag=plan->getMesh()->findTag("global_id");
  for (int dimension = 0; dimension <= plan->ghost_dim; ++dimension)
  {
    PCU_Comm_Begin();
    sendGhosts(plan, dimension, senders[dimension], tags);
    PCU_Comm_Send();
    EntityVector received;
    receiveGhosts(plan,global_id_tag,tags,received);
    setupGhosts(plan->getMesh(),received,senders[dimension]);
  }
}

// *********************************************************
void pumi_ghost_create(pMesh m, Ghosting* plan)
// *********************************************************
{
  if (PCU_Comm_Peers()==1) return;
  
  EntityVector affected[4];
  ghost_getAffected(m,plan,affected);
  EntityVector senders[4];
  getSenders(m,affected,senders);
//  reduceMatchingToSenders(m,senders);
//  updateResidences(m,plan,affected);
  ghost_unifyPids(plan, affected);
  ghost_moveEntities(plan, senders);
  delete plan;
//  updateMatching(m,affected,senders);
//  deleteOldEntities(m,affected);
  m->acceptChanges();
  // update global id

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
  pTag parts_index_tag=m->findTag("_parts_index_");
  if (parts_index_tag) m->destroyTag(parts_index_tag);
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

