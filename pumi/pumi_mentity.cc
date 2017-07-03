/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>

using std::vector;
// mesh entity functions
int pumi_ment_getDim(pMeshEnt e)
{
  return apf::getDimension(pumi::instance()->mesh, e);
}

int pumi_ment_getTopo(pMeshEnt e) 
{
  return pumi::instance()->mesh->getType(e);
}

int pumi_ment_getNumAdj(pMeshEnt e, int target_dim)
{
  int ent_dim= apf::getDimension(pumi::instance()->mesh, e);
  if (ent_dim<target_dim) // upward
  {
    apf::Adjacent adjacent;
    pumi::instance()->mesh->getAdjacent(e,target_dim,adjacent);      
    return adjacent.getSize();
  }
  else if (ent_dim>target_dim)
  {
    apf::Downward adjacent;
    return pumi::instance()->mesh->getDownward(e,target_dim,adjacent); 
  }
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<": invalid target dimension "<<target_dim<<"\n";
  return 0;
}

// if target_dim=-1, get all downward adjacent entities
void pumi_ment_getAdj(pMeshEnt e, int target_dim, std::vector<pMeshEnt>& vecAdjEnt)
{
  int ent_dim= apf::getDimension(pumi::instance()->mesh, e);
  if (ent_dim<target_dim) // upward
  {
    apf::Adjacent adjacent;
    pumi::instance()->mesh->getAdjacent(e,target_dim,adjacent);      
    for (size_t i=0; i<adjacent.getSize(); ++i)
      vecAdjEnt.push_back(adjacent[i]);
  }
  else if (ent_dim>target_dim)
  {
    if (target_dim>=0)
    {
      apf::Downward adjacent;
      int num_adj= pumi::instance()->mesh->getDownward(e,target_dim,adjacent); 
      for (int j = 0; j < num_adj; ++j) 
        vecAdjEnt.push_back(adjacent[j]);
    }
    else // get all downward adjacent entities
    {
      switch (ent_dim)
      {
        case 3:
          {
            apf::Downward faces;
            int num_faces= pumi::instance()->mesh->getDownward(e,2,faces); 
            for (int i=0; i<num_faces;++i)
            {
              vecAdjEnt.push_back(faces[i]);
              apf::Downward edges;
              int num_edges= pumi::instance()->mesh->getDownward(faces[i],1,edges); 
              for (int j=0; j<num_edges;++j)
                if (std::find(vecAdjEnt.begin(), vecAdjEnt.end(), edges[j])==vecAdjEnt.end())
                  vecAdjEnt.push_back(edges[j]);
            }
            apf::Downward vertices;
            int num_vertices= pumi::instance()->mesh->getDownward(e,0,vertices);
            for (int j=0; j<num_vertices;++j)
              vecAdjEnt.push_back(vertices[j]);
            if (num_vertices==4) assert(vecAdjEnt.size()==14);
            break;
          }
        case 2:// face
          {
            apf::Downward edges;
            int num_edges= pumi::instance()->mesh->getDownward(e,1,edges); 
            for (int j=0; j<num_edges;++j)
              vecAdjEnt.push_back(edges[j]);
            apf::Downward vertices;
            int num_vertices= pumi::instance()->mesh->getDownward(e,0,vertices);
            for (int j=0; j<num_vertices;++j)
              vecAdjEnt.push_back(vertices[j]);
            break;
          }
        case 1: // edge
          {
            apf::Downward vertices;
            int num_vertices= pumi::instance()->mesh->getDownward(e,0,vertices);
            for (int j=0; j<num_vertices;++j)
              vecAdjEnt.push_back(vertices[j]);
            break;
          }
        default: break;
      } // switch
    } // else
  } // else if
}

void pumi_ment_get2ndAdj (pMeshEnt e, int bridge_dim, int target_dim, std::vector<pMeshEnt>& vecAdjEnt)
{
  if (bridge_dim==target_dim) 
  {
    if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<": invalid bridge/target dimension \n";
    return;  // error
  }

  apf::Adjacent adjacent;
  apf::getBridgeAdjacent(pumi::instance()->mesh, e, bridge_dim, target_dim, adjacent);
  for (size_t i=0; i<adjacent.getSize(); ++i)
    vecAdjEnt.push_back(adjacent[i]);
}

int pumi_ment_getID(pMeshEnt e)
{
  return getMdsIndex(pumi::instance()->mesh, e);
}

pGeomEnt pumi_ment_getGeomClas(pMeshEnt e)
{
  gmi_ent* clas=(gmi_ent*)pumi::instance()->mesh->toModel(e);
  int dim = gmi_dim(pumi::instance()->model->getGmi(), clas);
  return pumi::instance()->model->getGeomEnt(dim, clas);
}

pMeshEnt pumi_medge_getOtherVtx(pMeshEnt edge, pMeshEnt vtx)
{
  return apf::getEdgeVertOppositeVert(pumi::instance()->mesh, edge, vtx);
}

// owner part information
// FIXME: based on sharing
int pumi_ment_getOwnPID(pMeshEnt e, pSharing shr)
{
  if (!shr)
    return pumi::instance()->mesh->getOwner(e);
  return shr->getOwner(e);
  
}

// FIXME: based on sharing
pMeshEnt pumi_ment_getOwnEnt(pMeshEnt e, pSharing shr)
{
  if (!(pumi::instance()->mesh->isShared(e))) // internal ent
    return e;
  int own_partid;
  if (!shr)
    own_partid= pumi::instance()->mesh->getOwner(e);
  else
    own_partid=shr->getOwner(e);

  if (own_partid==pumi_rank()) return e;
  return pumi_ment_getRmt(e, own_partid);
}

bool pumi_ment_isOwned(pMeshEnt e, pSharing shr)
{  
  if (!shr) 
    return (pumi_ment_getOwnPID(e)==pumi_rank());
  return shr->isOwned(e);
}

bool pumi_ment_isOn(pMeshEnt e, int partID)
{
  // FIXME: include ghost copy
  if (partID==pumi_rank()) return true;
  apf::Copies remotes;
  pumi::instance()->mesh->getRemotes(e,remotes);
  APF_ITERATE(Copies,remotes,rit)
    if (rit->first==partID) return true;
  apf::Copies ghosts;
  pumi::instance()->mesh->getGhosts(e,ghosts);
  APF_ITERATE(Copies,ghosts,rit)
    if (rit->first==partID) return true;
 return false;
}

// remote copy information
bool pumi_ment_isOnBdry (pMeshEnt e)
{
  if (pumi::instance()->mesh->isShared(e)) // internal ent
    return true;
  else 
    return false;
}

int pumi_ment_getNumRmt (pMeshEnt e)
{
  if (!(pumi::instance()->mesh->isShared(e))) // internal ent
    return 0;
  apf::Copies remotes;
  pumi::instance()->mesh->getRemotes(e,remotes);
  return remotes.size();
}

void pumi_ment_getAllRmt(pMeshEnt e, Copies& remotes)
{
  if (pumi::instance()->mesh->isShared(e))
    pumi::instance()->mesh->getRemotes(e,remotes);
}

pMeshEnt pumi_ment_getRmt(pMeshEnt& e, int pid)
{
  if (!(pumi::instance()->mesh->isShared(e))) // internal ent
    return NULL;

  Copies remotes;
  pumi::instance()->mesh->getRemotes(e,remotes);
  return remotes[pid];
}

void pumi_ment_setRmt(pMeshEnt, int, pMeshEnt)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

void pumi_ment_deleteRmt (pMeshEnt, int)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

void pumi_ment_cleanRmt (pMeshEnt)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

void pumi_ment_setPtnTopology (pMeshEnt)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

int pumi_ment_getGlobalID(pMeshEnt e)
{
  pMeshTag tag = pumi::instance()->mesh->findTag("global_id");
  if (!tag) return -1;
  int global_id;
  pumi::instance()->mesh->getIntTag(e, tag, &global_id);
  return global_id;
  
}

void pumi_ment_getResidence(pMeshEnt e, Parts& residence)
{
  // apf::getResidence takes part IDs from partiton classification
  // note ghosting part ID's are not stored apf::getResidence
  pumi::instance()->mesh->getResidence(e,residence);
  // now add ghost part id's 
  apf::Copies ghosts;
  pumi::instance()->mesh->getGhosts(e,ghosts);
  APF_ITERATE(Copies,ghosts,rit)
   residence.insert(rit->first);
}

void pumi_ment_getClosureResidence(pMeshEnt e, Parts& residence)
{
  pumi::instance()->mesh->getResidence(e,residence);

  apf::Downward vertices;
  int nvtx = pumi::instance()->mesh->getDownward(e,0,vertices);
  pMeshEnt vtx;
  
  for (int i=0; i<nvtx; ++i)
  {
    vtx=vertices[i];
    Parts vtx_res;
    pumi::instance()->mesh->getResidence(vtx,vtx_res);
    apf::unite(residence, vtx_res);
  }
}

// ghosting information
bool pumi_ment_isGhost(pMeshEnt e)
{
  return pumi::instance()->mesh->isGhost(e);
}

bool pumi_ment_isGhosted (pMeshEnt e)
{
  return pumi::instance()->mesh->isGhosted(e);
}

int pumi_ment_getNumGhost (pMeshEnt e)
{
  if (!pumi::instance()->mesh->isGhosted(e)) return 0;
  apf::Copies ghosts;
  return pumi::instance()->mesh->getGhosts(e,ghosts);
}

void pumi_ment_getAllGhost(pMeshEnt e, Copies& ghosts)
{
  if (pumi::instance()->mesh->isGhosted(e))
    pumi::instance()->mesh->getGhosts(e,ghosts);
}

pMeshEnt pumi_ment_getGhost(pMeshEnt& e, int pid)
{
 if (!(pumi::instance()->mesh->isGhosted(e))) // internal ent
    return NULL;

  Copies ghosts;
  pumi::instance()->mesh->getGhosts(e,ghosts);
  return ghosts[pid];
}

void pumi_ment_setNumber(pMeshEnt e, pNumbering n, int node, int component, int number)
{
  apf::number(n, e, node, component, number);
}

int pumi_ment_getNumber(pMeshEnt e, pNumbering n, int node, int component)
{
  return apf::getNumber(n, e, node, component);
}

bool pumi_ment_isNumbered(pMeshEnt e, pNumbering n)
{
  return apf::isNumbered(n, e, 0, 0);
}
