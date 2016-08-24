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
  assert (parts_index_tag);
  if (!PCU_Comm_Self()) std::cout<<__func__<<": deleting parts_vec_index from entities\n";
  // FIXME: this is not efficient
  for (int d=0; d<=ghost_dim; ++d)
    apf::removeTagFromDimension(m, parts_index_tag, d);
  m->destroyTag(parts_index_tag);
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

#include <assert.h>
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

  /* store the sender as a ghost copy */
  plan->getMesh()->addGhost(entity, from, sender);
//  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<entity<<"(d "<<getDimension(plan->getMesh(), entity)
//           <<", id "<<global_id<<")->addGhost("<<from<<", "<<sender<<")\n";
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
    if (!plan->has(e))
      continue;

    Copies remotes;
    plan->getMesh()->getRemotes(e,remotes);

    Parts sendTo;
    apf::split(remotes,plan->sending(e, ent_dim),sendTo);
    if (plan->getMesh()->isGhosted(e)) 
    {
      Copies ghosts;
      plan->getMesh()->getGhosts(e,ghosts);
      apf::split(ghosts,plan->sending(e, ent_dim),sendTo);
    }

    if (!sendTo.size())
      continue;
    APF_ITERATE(Parts,sendTo,sit)
      apf::packEntity(plan->getMesh(),*sit,e,tags);
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
//    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<entity<<"(d "<<getDimension(m,entity)
//             <<", id "<<pumi_ment_getglobalid(entity)<<")->addGhost("<<from<<", "<<sender<<")\n";
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

// *****************************************
static void ghost_collectEntities (pMesh m, Ghosting* plan, EntityVector entitiesToGhost[4])
// *****************************************
{    

  pMeshEnt down_ent; 
  pMeshEnt remote_ent;
  pMeshEnt ghost_ent;
  int dummy=1;
  std::vector<pMeshEnt> DownEnts;
  DownEnts.resize(27);
  
  int down_ent_dim, ghost_dim = plan->ghost_dim;

  pTag tag = m->findTag("entity_2_ghost");
  if (!tag)
    tag = m->createIntTag("entity_2_ghost",1);

  entitiesToGhost[ghost_dim].reserve(plan->count());

  pMeshEnt e;
  apf::MeshIterator* it = m->begin(ghost_dim);
  while ((e = m->iterate(it)))
  {
    if (plan->has(e))
      entitiesToGhost[ghost_dim].push_back(e);
  }
  m->end(it);

  // SEOL: this should be fixed
  APF_ITERATE(EntityVector,entitiesToGhost[ghost_dim],it)
  {
    ghost_ent=*it;    
//    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": -- ghost_ent "<<pumi_ment_getglobalid(ghost_ent)<<"\n";
    DownEnts.clear();
    pumi_ment_getadj(ghost_ent, -1, DownEnts);
  
    for (std::vector<pMeshEnt>::iterator downadj_it=DownEnts.begin();downadj_it!=DownEnts.end();++downadj_it)
    {
      down_ent = *downadj_it;
      down_ent_dim=getDimension(m, down_ent);
      if (!m->hasTag(down_ent,tag))
      {
        m->setIntTag(down_ent,tag,&dummy);
        entitiesToGhost[down_ent_dim].push_back(down_ent);
//        std::cout<<"("<<pumi_rank()<<") "<<__func__<<": entitiesToGhost["<<down_ent_dim<<"].push_back("<<pumi_ment_getglobalid(down_ent)<<")\n";
      }
      APF_ITERATE(Parts, plan->sending(ghost_ent, ghost_dim), pit)
        plan->send(down_ent, *pit);
    } // for (std::vector<pMeshEnt>::iterator downadj_it
  } // APF_ITERATE

  // do communication to unify ghost target pids
  void* msg_send;
  pMeshEnt* s_ent;
  int* s_id; 

  size_t msg_size;
  int numBP;
  for (int dim = 0; dim <=ghost_dim; ++dim)
  {
    PCU_Comm_Begin();
    APF_ITERATE(EntityVector,entitiesToGhost[dim],it)
    {
      e = *it;
      if (!m->isShared(e)) continue;

      int num_pids=plan->count(e, dim);
      apf::Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(apf::Copies,remotes,rit)
      {
        msg_size=sizeof(pMeshEnt) +num_pids*sizeof(int);
        msg_send = malloc(msg_size);
      
        s_ent = (pMeshEnt*)msg_send; 
        *s_ent = rit->second; 
        int *pids = (int*)((char*)msg_send + sizeof(pMeshEnt));
        int pos = 0;
       
        APF_ITERATE(Parts, plan->sending(e, dim), pit)
        {
          pids[pos]=*pit;
          ++pos;
        }
        PCU_Comm_Write(rit->first, (void*)msg_send, msg_size);
        free(msg_send);    
      }
    } // for entitiesToGhost[dim]
    PCU_Comm_Send();
  
    // receive phase
    void *msg_recv;
    int pid_from, r_dim, r_id;
    int* pids;
    pMeshEnt r;
    while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
    {
      r = *((pMeshEnt*)msg_recv); 
      if ( !m->hasTag(r,tag))
      {
        m->setIntTag(r,tag,&dummy);
        entitiesToGhost[dim].push_back(r);
//        std::cout<<"("<<pumi_rank()<<") "<<__func__<<": rmt entitiesToGhost["<<dim<<"].push_back("<<pumi_ment_getglobalid(r)<<")\n";
      }

      pids = (int*)((char*)msg_recv+sizeof(pMeshEnt)); 
      int num_pids = (msg_size-sizeof(pMeshEnt))/sizeof(int);
      for (int i = 0; i < num_pids; ++i)
        plan->send(r, pids[i]); //parts_vec[r_dim][r_id].insert(pids[i]);
    } // while

    APF_ITERATE(EntityVector,entitiesToGhost[dim],it)
      m->removeTag(*it,tag);
  } // for dim

  m->destroyTag(tag);
}
// *********************************************************
void pumi_ghost_create(pMesh m, Ghosting* plan)
// *********************************************************
{
  if (PCU_Comm_Peers()==1) return;
  
  EntityVector affected[4];
  ghost_collectEntities(m, plan, affected);
  EntityVector senders[4];
  getSenders(m,affected,senders);
  ghost_moveEntities(plan, senders);

  delete plan;
  m->acceptChanges();
}


// *********************************************************
void pumi_ghost_createlayer (pMesh m, int brg_dim, int ghost_dim, int num_layer, int include_copy)
// *********************************************************
{
  if (PCU_Comm_Peers()==1 || num_layer==0) return;
  
  int dummy=1, mesh_dim=m->getDimension(), self = pumi_rank();;
  
// brid/ghost dim check
  if (brg_dim>=ghost_dim || 0>brg_dim || brg_dim>=mesh_dim || 
      ghost_dim>mesh_dim || ghost_dim<1)
  {
    if (!self)
       std::cout<<__func__<<" ERROR: invalid bridge/ghost dimension\n";   
    return;
  }
  pTag tag = m->createIntTag("entity_2_ghost",1);
  Ghosting* plan = new Ghosting(m, ghost_dim);

// ********************************************
// STEP 1: compute entities to ghost
// ********************************************

  pMeshEnt ghost_ent;
  pMeshEnt brg_ent;

  std::vector<pMeshEnt> processed_ent;
  std::vector<pMeshEnt> adj_ent;
  //PUMI_PartEntIter_InitPartBdry(part, PUMI_ALL, brg_dim, PUMI_ALLTOPO, pbdry_iter);
  //while (PUMI_PartEntIter_GetNext(pbdry_iter, brg_ent)==PUMI_SUCCESS)

  apf::MeshIterator* it = m->begin(brg_dim);
  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<__LINE__<<"\n";

  while ((brg_ent = m->iterate(it)))
  {
    if (!m->isShared(brg_ent)) continue; // skip non-partboundary entity
    if (!include_copy && m->getOwner(brg_ent)!=self) continue;

    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<__LINE__<<" brg_ent "<<pumi_ment_getglobalid(brg_ent)<<"\n";
    processed_ent.clear();

    apf::Adjacent adjacent;
    m->getAdjacent(brg_ent,ghost_dim, adjacent);   
    APF_ITERATE(apf::Adjacent, adjacent, adj_ent_it)
    {
      ghost_ent = *adj_ent_it;
      if (m->isGhost(ghost_ent)) continue; // skip ghost copy

      //copy_RC_to_BP(brg_ent, ghost_ent);
      apf::Copies remotes;
      m->getRemotes(brg_ent,remotes);
      APF_ITERATE(apf::Copies,remotes,rit)
        plan->send(ghost_ent, rit->first);

      m->setIntTag(ghost_ent,tag,&dummy);
      processed_ent.push_back(ghost_ent);
    }
    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<__LINE__<<"\n";
    
    //for (std::vector<pMeshEnt>::iterator ghost_it=brg_ghost_map[ent]->begin();ghost_it!=brg_ghost_map[ent]->end(); ++ghost_it)
//    std::cout<<"[p"<<PUMI_CommRank()<<"] "<<PUMI_MeshEnt_StrID(brg_ent)<<"- 1-layer ghost "<<processed_ent.size()<<"\n";
    int start_prev_layer=0, size_prev_layer=processed_ent.size(), num_prev_layer;
    for (int layer=2; layer<num_layer+1; ++layer)
    {  
      num_prev_layer=0;
      for (int i=start_prev_layer; i<size_prev_layer; ++i)
      {
        ghost_ent = processed_ent.at(i);
        adj_ent.clear();
        pumi_ment_get2ndadj (ghost_ent, brg_dim, ghost_dim, adj_ent);
        
        for (std::vector<pMeshEnt>::iterator adj_ent_it=adj_ent.begin(); adj_ent_it!=adj_ent.end(); ++adj_ent_it)
        {
          if (m->isGhost(*adj_ent_it) || m->hasTag(*adj_ent_it,tag))
            continue; // skip ghost copy or already-processed copy
      
          //copy_RC_to_BP(brg_ent, ghost_ent);
          apf::Copies remotes;
          m->getRemotes(brg_ent,remotes);
          APF_ITERATE(apf::Copies,remotes,rit)
            plan->send(*adj_ent_it, rit->first);

          m->setIntTag(*adj_ent_it,tag,&dummy);
          processed_ent.push_back(*adj_ent_it);
          ++num_prev_layer;
        } // for (std::vector<pMeshEnt>::iterator adj_ent_it=adj_ent.begin()
      } // for int i=start_prev_layerghost_it
      start_prev_layer+=size_prev_layer;
      size_prev_layer+=num_prev_layer;
    } // for layer
    std::cout<<"("<<pumi_rank()<<") "<<__func__<<": "<<__LINE__<<" brg_ent "<<pumi_ment_getglobalid(brg_ent)<<"\n";
  } // while brg_ent
  m->end(it);

// ********************************************
// STEP 2: perform ghosting
// ********************************************
  std::cout<<"("<<pumi_rank()<<") "<<__func__<<": plan->count()="<<plan->count()<<"\n";

  pumi_ghost_create(m, plan);
}

// *********************************************************
void pumi_ghost_delete (pMesh m)
// *********************************************************
{
  pMeshEnt e;
  for (int d=4; d>=0; --d)
  {
    apf::MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it)))
    {
      if (m->isGhosted(e))
        m->deleteGhost(e);
      if (m->isGhost(e))
        m->destroy(e);
    }
    m->end(it);
  }
}

// *********************************************************
void pumi_ghost_info (pMesh m, std::vector<int>& ghostinfo)
// *********************************************************
{
  if (!pumi_rank()) 
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: not supported\n";
}

