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
#include <pcu_util.h>
#include <cstdlib>

#include "apf.h"
#include "apfMDS.h"

using std::map;
using std::set;

Ghosting::Ghosting(pMesh mesh, int d)
{
  m = mesh;
  ghost_dim = d;

  if (!m->findTag("ghost_tag"))
    pumi::instance()->ghost_tag = m->createIntTag("ghost_tag",1);
  if (!m->findTag("ghosted_tag"))
    pumi::instance()->ghosted_tag = m->createIntTag("ghosted_tag",1);

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
  PCU_ALWAYS_ASSERT(parts_index_tag);

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
//  if (to==PCU_Comm_Self()) return;

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
  PCU_ALWAYS_ASSERT(index!=-1);
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
      std::cout<<"("<<PCU_Comm_Self()<<") ghost e "<<getMdsIndex(m,e)<<" to "<<*pit<<"\n";
  }
  m->end(it);
}

Parts& Ghosting::sending(pMeshEnt e, int d)
{
  int index;
  PCU_ALWAYS_ASSERT(m->hasTag(e, parts_index_tag));
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

#include <pcu_util.h>
// *********************************************************
static pMeshEnt unpackGhost(Ghosting* plan, apf::DynamicArray<pMeshTag>& tags)
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
  residence.insert(from);
  plan->getMesh()->setResidence(entity,residence);
  apf::unpackTags(plan->getMesh(),entity,tags);
  apf::unpackRemotes(plan->getMesh(),entity);

  /* store the sender as a ghost copy */
  plan->getMesh()->addGhost(entity, from, sender);
  pumi::instance()->ghost_vec[apf::getDimension(plan->getMesh(), entity)].push_back(entity);
  plan->getMesh()->setIntTag(entity,pumi::instance()->ghost_tag,&from);
  return entity;
}

// *********************************************************
static void ghost_receiveEntities(Ghosting* plan, apf::DynamicArray<pMeshTag>& tags,
    EntityVector& received)
// *********************************************************
{
  received.reserve(1024);
  while (PCU_Comm_Receive())
    received.push_back(unpackGhost(plan,tags));
}

// *********************************************************
static void setupGhosts(pMesh m, EntityVector& received)
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
    apf::Copies remotes;
    m->getRemotes(entity, remotes);
    APF_ITERATE(Copies,remotes,rit)
    {
      PCU_COMM_PACK(rit->first, rit->second); // sender
      PCU_COMM_PACK(rit->first,entity);
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
  {
    int from = PCU_Comm_Sender();
    pMeshEnt entity;
    PCU_COMM_UNPACK(entity);
    pMeshEnt sender;
    PCU_COMM_UNPACK(sender);
    m->addGhost(entity, from, sender);
    if (!m->hasTag(entity, pumi::instance()->ghosted_tag))
    {
      pumi::instance()->ghosted_vec[apf::getDimension(m, entity)].push_back(entity);
      m->setIntTag(entity,pumi::instance()->ghosted_tag,&from);
    }
  }
}

// *****************************************
static void ghost_collectEntities (pMesh m, Ghosting* plan, EntityVector entitiesToGhost[4])
// *****************************************
{    

  pMeshEnt down_ent; 
  pMeshEnt ghost_ent;
  int dummy=1;
  std::vector<pMeshEnt> DownEnts;
  DownEnts.resize(27);
  
  int down_ent_dim, ghost_dim = plan->ghost_dim;

  pMeshTag tag = m->findTag("entity_2_ghost");
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
    DownEnts.clear();
    pumi_ment_getAdj(ghost_ent, -1, DownEnts);
  
    for (std::vector<pMeshEnt>::iterator downadj_it=DownEnts.begin();downadj_it!=DownEnts.end();++downadj_it)
    {
      down_ent = *downadj_it;
      down_ent_dim=getDimension(m, down_ent);
      if (!m->hasTag(down_ent,tag))
      {
        m->setIntTag(down_ent,tag,&dummy);
        entitiesToGhost[down_ent_dim].push_back(down_ent);
      }
      APF_ITERATE(Parts, plan->sending(ghost_ent, ghost_dim), pit)
        plan->send(down_ent, *pit);
    } // for (std::vector<pMeshEnt>::iterator downadj_it
  } // APF_ITERATE

  // do communication to unify ghost target pids
  void* msg_send;
  pMeshEnt* s_ent;

  size_t msg_size;
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
    int pid_from;
    int* pids;
    pMeshEnt r;
    while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
    {
      r = *((pMeshEnt*)msg_recv); 
      if ( !m->hasTag(r,tag))
      {
        m->setIntTag(r,tag,&dummy);
        entitiesToGhost[dim].push_back(r);
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

// *****************************************
void set_subtract(std::set<int> A, std::set<int> B, std::set<int>& C)
// *****************************************
{
  std::set<int>::iterator aiter, biter;
  for (aiter=A.begin(); aiter!=A.end();++aiter)
  {
    biter = B.find(*aiter);
    if (biter==B.end()) C.insert(*aiter);
  }
} 


// **********************************************
void ghost_sendEntities(Ghosting* plan, int entDim,
      std::vector<pMeshEnt>& entitiesToExchg, apf::DynamicArray<pMeshTag>& tags)	  
// **********************************************
{
  pMeshEnt ent;
  int src_partid=PCU_Comm_Self();
  pMesh m = plan->getMesh();

  std::set<int> res_parts, temp; 
  std::set<int> target_pids;
  for(std::vector<pMeshEnt>::const_iterator eit=entitiesToExchg.begin(); eit!=entitiesToExchg.end();++eit)
  {
    ent = *eit;

    if (src_partid!=m->getOwner(ent)) continue;

    res_parts.clear(); 
    temp.clear();
    target_pids.clear();
 
    res_parts.insert(src_partid);

    if (m->isShared(ent)) // let the owner part send the ghost copy
    {
      apf::Copies remotes;
      m->getRemotes(ent,remotes);
      APF_ITERATE(apf::Copies,remotes,rit)
        res_parts.insert(rit->first);
    }

    if (m->isGhosted(ent))
    {
      apf::Copies ghosts;
      m->getGhosts(ent,ghosts);
      APF_ITERATE(apf::Copies,ghosts,rit)
        res_parts.insert(rit->first);
    }

    APF_ITERATE(Parts, plan->sending(ent, entDim), pit)
      target_pids.insert(*pit);
    
    set_subtract(target_pids,res_parts, temp);
    if (temp.size()==0) continue;

    for (std::set<int>::iterator piter=temp.begin(); piter!=temp.end();++piter)
    {
      if (*piter==src_partid) continue;
      apf::packEntity(plan->getMesh(),*piter,ent,tags, true/*ghosting*/);  
    }
  }
}

#include "apfNumbering.h"
#include "apfShape.h"
// *********************************************************
void pumi_ghost_create(pMesh m, Ghosting* plan)
// *********************************************************
{
  if (PCU_Comm_Peers()==1) return;
 
  std::vector<apf::Field*> fields;
  std::vector<apf::Field*> frozen_fields;
  for (int i=0; i<m->countFields(); ++i)
  {
    pField f = m->getField(i);
    if (isFrozen(f))
    {
      frozen_fields.push_back(f); // turn field data from tag to array
      apf::unfreeze(f);
    }
  }

  double t0=PCU_Time();

  EntityVector entities_to_ghost[4];
  ghost_collectEntities(m, plan, entities_to_ghost);

  apf::DynamicArray<pMeshTag> tags;
  plan->getMesh()->getTags(tags);
  for (int dimension = 0; dimension <= plan->ghost_dim; ++dimension)
  {
    PCU_Comm_Begin();
    ghost_sendEntities(plan, dimension, entities_to_ghost[dimension], tags);
    PCU_Comm_Send();
    EntityVector received;
    ghost_receiveEntities(plan,tags,received);
    setupGhosts(plan->getMesh(),received);
  }
  
  delete plan;
  m->acceptChanges();
  
  while (m->countNumberings()) destroyNumbering(m->getNumbering(0));
  for (std::vector<apf::Field*>::iterator fit=frozen_fields.begin(); fit!=frozen_fields.end(); ++fit)
    apf::freeze(*fit);    

  if (!PCU_Comm_Self())
    printf("mesh ghosted in %f seconds\n", PCU_Time()-t0);
}

// *********************************************************
void do_off_part_bridge(pMesh m, int brg_dim, int ghost_dim, int num_layer, 
                        std::set<pMeshEnt>** off_bridge_set, Ghosting* plan)
// *********************************************************
{
  std::map<pMeshEnt, set<int> > off_bridge_marker;
  pMeshEnt ghost_ent;
  pMeshEnt brg_ent;
  void* msg_send;
  pMeshEnt* s_ent;
  size_t msg_size;
  int dummy=1;
  PCU_Comm_Begin();

  for (int layer=2; layer<num_layer+1; ++layer)
  {
    for (int pid=0; pid<pumi_size(); ++pid)
    {
      for (std::set<pMeshEnt>::iterator off_it=off_bridge_set[layer][pid].begin(); 
            off_it!=off_bridge_set[layer][pid].end(); ++off_it)
      {
        brg_ent = *off_it;
        apf::Copies brg_remotes;
        m->getRemotes(brg_ent,brg_remotes);
        APF_ITERATE(apf::Copies,brg_remotes,brg_rit)
        {
          if (brg_rit->first==pid) continue;
          msg_size=sizeof(pMeshEnt) +2*sizeof(int);
          msg_send = malloc(msg_size);
          s_ent = (pMeshEnt*)msg_send; 
          *s_ent = brg_rit->second; 
          int *s_int = (int*)((char*)msg_send + sizeof(pMeshEnt));
          s_int[0]=layer;
          s_int[1]=pid;
          PCU_Comm_Write(brg_rit->first, (void*)msg_send, msg_size);
          free(msg_send);    
        }  // APF_ITERATE
      } // for (std::map<pMeshEnt, set<int> >::iterator iter     
    }
  }
  // clean up
  for (int i=0; i<num_layer+1;++i)
    for (int j=0; j<pumi_size();++j)
      off_bridge_set[i][j].clear();

  PCU_Comm_Send();
  // receive phase
  void *msg_recv;
  int pid_from;
  int* r_int;
  pMeshEnt r;
  int r_layer, r_pid;
  pMeshTag tag = m->findTag("ghost_check_mark");
  std::vector<pMeshEnt> processed_ent;
  std::vector<pMeshEnt> adj_ent;

  while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
  {
    processed_ent.clear();

    r = *((pMeshEnt*)msg_recv); 
    r_int = (int*)((char*)msg_recv+sizeof(pMeshEnt)); 
    r_layer=r_int[0];
    r_pid=r_int[1];
    off_bridge_marker[r].insert(r_pid);
    apf::Adjacent ghost_cands;
    m->getAdjacent(r,ghost_dim, ghost_cands);   
    APF_ITERATE(apf::Adjacent, ghost_cands, adj_ent_it)
    {
      ghost_ent = *adj_ent_it;
      if (m->isGhost(ghost_ent)) continue; // skip ghost copy
      plan->send(ghost_ent, r_pid);

      m->setIntTag(ghost_ent,tag,&dummy);
      processed_ent.push_back(ghost_ent);
    
      if (r_layer<num_layer)
      {
        apf::Downward adjacent;
        int num_brg=m->getDownward(ghost_ent,brg_dim, adjacent);
        for (int b=0; b<num_brg; ++b)
        {     
          if (m->isShared(adjacent[b]) && adjacent[b]!=r && !pumi_ment_isOn(adjacent[b], r_pid))
            off_bridge_set[r_layer+1][r_pid].insert(adjacent[b]);
        }
      } // if (r_layer<num_layer)
    } // APF_ITERATE

    int start_prev_layer=0, size_prev_layer=processed_ent.size(), num_prev_layer;
    for (int layer=r_layer+1; layer<num_layer+1; ++layer)
    {  
      num_prev_layer=0;
      for (int i=start_prev_layer; i<size_prev_layer; ++i)
      {
        ghost_ent = processed_ent.at(i);
        adj_ent.clear();
        pumi_ment_get2ndAdj (ghost_ent, brg_dim, ghost_dim, adj_ent);

        for (std::vector<pMeshEnt>::iterator git=adj_ent.begin(); git!=adj_ent.end(); ++git)
        {
          if (m->isGhost(*git) || m->hasTag(*git,tag))
            continue; // skip ghost copy or already-processed copy
      
          plan->send(*git, r_pid);

          m->setIntTag(*git,tag,&dummy);
          processed_ent.push_back(*git);
          ++num_prev_layer;
        } // for (std::vector<pMeshEnt>::iterator git=adj_ent.begin()

        // collect off-part adjacent bridges
        if (layer<num_layer)
        {
          apf::Downward adjacent;
          int num_brg=m->getDownward(ghost_ent,brg_dim, adjacent);
          for (int b=0; b<num_brg; ++b)
          {     
            if (m->isShared(adjacent[b]) && adjacent[b]!=r && !pumi_ment_isOn(adjacent[b], r_pid))
              off_bridge_set[layer+1][r_pid].insert(adjacent[b]);
          } // for int b=0
        } // if (layer<=num_layer)
      } // for int i=start_prev_layer
      start_prev_layer+=size_prev_layer;
      size_prev_layer+=num_prev_layer;
    } // for layer
    for (std::vector<pMeshEnt>::iterator git=processed_ent.begin(); git!=processed_ent.end(); ++git)
      m->removeTag(*git,tag);
  } // while r

  if (num_layer==2) return;

  // the below runs only if num_layer>2
  int global_num_off_part, local_num_off_part=0;
  for (int i=0; i<num_layer+1;++i)
    for (int j=0; j<pumi_size();++j)
      local_num_off_part+=off_bridge_set[i][j].size();
  MPI_Allreduce(&local_num_off_part, &global_num_off_part, 1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  while (global_num_off_part)
  {
    PCU_Comm_Begin();

    for (int layer=0; layer<num_layer+1; ++layer)
    {
      for (int pid=0; pid<pumi_size(); ++pid)
      {
        for (std::set<pMeshEnt>::iterator off_it=off_bridge_set[layer][pid].begin(); 
            off_it!=off_bridge_set[layer][pid].end(); ++off_it)
        {
          brg_ent = *off_it;
          apf::Copies brg_remotes;
          m->getRemotes(brg_ent,brg_remotes);
          APF_ITERATE(apf::Copies,brg_remotes,brg_rit)
          {
            if (brg_rit->first==pid) continue;
            msg_size=sizeof(pMeshEnt) +2*sizeof(int);
            msg_send = malloc(msg_size);
            s_ent = (pMeshEnt*)msg_send; 
            *s_ent = brg_rit->second; 
            int *p_int = (int*)((char*)msg_send + sizeof(pMeshEnt));
            p_int[0]=layer;
            p_int[1]=pid;
            PCU_Comm_Write(brg_rit->first, (void*)msg_send, msg_size);
            free(msg_send);    
          }  // APF_ITERATE
        } // for off_it     
      }
    }
    // clean up
    for (int i=0; i<num_layer+1;++i)
      for (int j=0; j<pumi_size();++j)
        off_bridge_set[i][j].clear();

    PCU_Comm_Send();

    while (PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
    {
      processed_ent.clear();

      r = *((pMeshEnt*)msg_recv); 
      r_int = (int*)((char*)msg_recv+sizeof(pMeshEnt)); 

      r_layer=r_int[0];
      r_pid=r_int[1];
      PCU_ALWAYS_ASSERT(r_layer<=num_layer && r_pid<pumi_size());
      off_bridge_marker[r].insert(r_pid);

      apf::Adjacent ghost_cands;
      m->getAdjacent(r,ghost_dim, ghost_cands);   
      APF_ITERATE(apf::Adjacent, ghost_cands, adj_ent_it)
      {
        ghost_ent = *adj_ent_it;
        if (m->isGhost(ghost_ent)) continue; // skip ghost copy
        plan->send(ghost_ent, r_pid);

        m->setIntTag(ghost_ent,tag,&dummy);
        processed_ent.push_back(ghost_ent);

        if (r_layer<num_layer)
        {
          apf::Downward adjacent;
          int num_brg=m->getDownward(ghost_ent,brg_dim, adjacent);
          for (int b=0; b<num_brg; ++b)
          {     
            if (m->isShared(adjacent[b]) && adjacent[b]!=r && !pumi_ment_isOn(adjacent[b], r_pid))
            {
              if (off_bridge_marker[adjacent[b]].find(r_pid)==off_bridge_marker[adjacent[b]].end())
                off_bridge_set[r_layer+1][r_pid].insert(adjacent[b]);
            }
          }
        } // if (layer+1<=num_layer)
      } // APF_ITERATE

      int start_prev_layer=0, size_prev_layer=processed_ent.size(), num_prev_layer;
      for (int layer=r_layer+1; layer<num_layer+1; ++layer)
      {  
        num_prev_layer=0;
        for (int i=start_prev_layer; i<size_prev_layer; ++i)
        {
          ghost_ent = processed_ent.at(i);
          adj_ent.clear();
          pumi_ment_get2ndAdj (ghost_ent, brg_dim, ghost_dim, adj_ent);
  
          for (std::vector<pMeshEnt>::iterator git=adj_ent.begin(); git!=adj_ent.end(); ++git)
          {
            if (m->isGhost(*git) || m->hasTag(*git,tag))
              continue; // skip ghost copy or already-processed copy
      
            plan->send(*git, r_pid);
 
            m->setIntTag(*git,tag,&dummy);
            processed_ent.push_back(*git);
            ++num_prev_layer;
          } // for (std::vector<pMeshEnt>::iterator git=adj_ent.begin()

          // collect off-part adjacent bridges
          if (layer<num_layer)
          {
            apf::Downward adjacent;
            int num_brg=m->getDownward(ghost_ent,brg_dim, adjacent);
            for (int b=0; b<num_brg; ++b)
            {     
              if (m->isShared(adjacent[b]) && adjacent[b]!=r && !pumi_ment_isOn(adjacent[b], r_pid))
                off_bridge_set[layer+1][r_pid].insert(adjacent[b]);
            } // for int b=0
          } // if (layer<=num_layer)
        } // for int i=start_prev_layer
        start_prev_layer+=size_prev_layer;
        size_prev_layer+=num_prev_layer;
      } // for layer
      for (std::vector<pMeshEnt>::iterator git=processed_ent.begin(); git!=processed_ent.end(); ++git)
        m->removeTag(*git,tag);
    }  // while (PCU_Comm_Read)
    local_num_off_part=0;
    for (int i=0; i<num_layer+1;++i)
      for (int j=0; j<pumi_size();++j)
        local_num_off_part+=off_bridge_set[i][j].size();
    MPI_Allreduce(&local_num_off_part, &global_num_off_part, 1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  } // while global_off_part_brg
}


// *********************************************************
void pumi_ghost_createLayer (pMesh m, int brg_dim, int ghost_dim, int num_layer, int include_copy)
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

  double t0 = PCU_Time();

  pMeshTag tag = m->createIntTag("ghost_check_mark",1);
  Ghosting* plan = new Ghosting(m, ghost_dim);

// ********************************************
// STEP 1: compute entities to ghost
// ********************************************

  pMeshEnt ghost_ent;
  pMeshEnt brg_ent;

  std::vector<pMeshEnt> processed_ent;
  std::vector<pMeshEnt> adj_ent;

  std::set<pMeshEnt>** off_bridge_set=new std::set<pMeshEnt>*[num_layer+1];
  for (int i=0; i<num_layer+1;++i)
    off_bridge_set[i]=new std::set<pMeshEnt>[pumi_size()];

  apf::MeshIterator* it = m->begin(brg_dim);
  while ((brg_ent = m->iterate(it)))
  {
    if (!m->isShared(brg_ent)) continue; // skip non-partboundary entity
    if (!include_copy && m->getOwner(brg_ent)!=self) continue;
    processed_ent.clear();

    apf::Copies remotes;
    m->getRemotes(brg_ent,remotes);

    apf::Adjacent adjacent;
    m->getAdjacent(brg_ent,ghost_dim, adjacent);   
    APF_ITERATE(apf::Adjacent, adjacent, adj_ent_it)
    {
      ghost_ent = *adj_ent_it;
      if (m->isGhost(ghost_ent)) continue; // skip ghost copy

      APF_ITERATE(apf::Copies,remotes,rit)
        plan->send(ghost_ent, rit->first);

      m->setIntTag(ghost_ent,tag,&dummy);
      processed_ent.push_back(ghost_ent);
    }

    int start_prev_layer=0, size_prev_layer=processed_ent.size(), num_prev_layer;
    for (int layer=2; layer<num_layer+1; ++layer)
    {  
      num_prev_layer=0;
      for (int i=start_prev_layer; i<size_prev_layer; ++i)
      {
        ghost_ent = processed_ent.at(i);
        adj_ent.clear();
        pumi_ment_get2ndAdj (ghost_ent, brg_dim, ghost_dim, adj_ent);

        for (std::vector<pMeshEnt>::iterator git=adj_ent.begin(); git!=adj_ent.end(); ++git)
        {
          if (m->isGhost(*git) || m->hasTag(*git,tag))
            continue; // skip ghost copy or already-processed copy
      
          APF_ITERATE(apf::Copies,remotes,rit)
            plan->send(*git, rit->first);

          m->setIntTag(*git,tag,&dummy);
          processed_ent.push_back(*git);
          ++num_prev_layer;
        } // for (std::vector<pMeshEnt>::iterator git=adj_ent.begin()

        // collect off-part adjacent bridges
        if (layer<=num_layer)
        {
          apf::Downward adjacent;
          int num_brg=m->getDownward(ghost_ent,brg_dim, adjacent);
          for (int b=0; b<num_brg; ++b)
          {     
            if (m->isShared(adjacent[b]) && brg_ent!=adjacent[b])
            {
              APF_ITERATE(apf::Copies,remotes,rit)
              {
                if (!pumi_ment_isOn(adjacent[b], rit->first))
                  off_bridge_set[layer][rit->first].insert(adjacent[b]);
              }
            } // if (m->isShared(adjacent[b])
          } // for int b=0
        } // if (layer<=num_layer)
      } // for int i=start_prev_layer

      start_prev_layer+=size_prev_layer;
      size_prev_layer+=num_prev_layer;
    } // for layer
    for (std::vector<pMeshEnt>::iterator git=processed_ent.begin(); git!=processed_ent.end(); ++git)
      m->removeTag(*git,tag);
  } // while brg_ent
  m->end(it);

// ********************************************
// STEP 2: deal with off-part adjacency, if any
// ********************************************
  int global_num_off_part=0, local_num_off_part=0;
  if (num_layer>=2)
  {
    for (int i=0; i<num_layer+1;++i)
      for (int j=0; j<pumi_size();++j)
        local_num_off_part+=off_bridge_set[i][j].size();
    MPI_Allreduce(&local_num_off_part, &global_num_off_part, 1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }

  if (global_num_off_part)
    do_off_part_bridge(m, brg_dim, ghost_dim, num_layer, off_bridge_set, plan);

  // clean up
  m->destroyTag(tag);
  for (int i=0; i<num_layer+1;++i)
    for (int j=0; j<pumi_size();++j)
      off_bridge_set[i][j].clear();

  for (int i=0; i<num_layer+1;++i)
    delete [] off_bridge_set[i];
  delete [] off_bridge_set;

// ********************************************
// STEP 3: perform ghosting
// ********************************************
  if (!PCU_Comm_Self())
    printf("ghosting plan computed in %f seconds\n", PCU_Time()-t0);

  pumi_ghost_create(m, plan);
}

// *********************************************************
void pumi_ghost_delete (pMesh m)
// *********************************************************
{
  pMeshTag tag = pumi::instance()->ghosted_tag;
  if (!tag) return;

  std::vector<apf::Field*> frozen_fields;
  for (int i=0; i<m->countFields(); ++i)
  {
    pField f = m->getField(i);
    if (isFrozen(f))
    {
      frozen_fields.push_back(f); // turn field data from tag to array
      apf::unfreeze(f);
    }
  }

  for (int d=3; d>=0; --d)
  {
    for (std::vector<pMeshEnt>::iterator it=pumi::instance()->ghost_vec[d].begin();
         it!=pumi::instance()->ghost_vec[d].end(); ++it)
      m->destroy(*it);
    for (std::vector<pMeshEnt>::iterator it=pumi::instance()->ghosted_vec[d].begin();
         it!=pumi::instance()->ghosted_vec[d].end(); ++it)
    {
      m->removeTag(*it,tag);
      m->deleteGhost(*it);
    }
  }

  // update owners
  m->acceptChanges();

  // delete tag
  m->destroyTag(pumi::instance()->ghost_tag);
  pumi::instance()->ghost_tag = NULL;
  m->destroyTag(pumi::instance()->ghosted_tag);
  pumi::instance()->ghosted_tag = NULL;

  for (int d=3; d>=0; --d)
  {
    pumi::instance()->ghost_vec[d].clear();
    pumi::instance()->ghosted_vec[d].clear();
  }

  while (m->countNumberings()) destroyNumbering(m->getNumbering(0));
  for (std::vector<apf::Field*>::iterator fit=frozen_fields.begin(); fit!=frozen_fields.end(); ++fit)
    apf::freeze(*fit);    
}

// *********************************************************
void pumi_ghost_getInfo (pMesh, std::vector<int>&)
// *********************************************************
{
  if (!pumi_rank()) 
    std::cout<<"[PUMI ERROR] "<<__func__<<" failed: not supported\n";
}

