/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <mpi.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>

using std::vector;
// geometric entity functions

int pumi_gent_getdim(pGeomEnt ge)
{
    return gmi_dim(pumi::instance()->model, ge);
}

int pumi_gent_getid(pGeomEnt ge)
{
  return gmi_tag(pumi::instance()->model, ge);
}

void pumi_gent_getrevclas (pGeomEnt g, std::vector<pMeshEnt>& ents)
{
  assert(!ents.size());
  int dim=gmi_dim(pumi::instance()->model, g);
  pMesh m = pumi::instance()->mesh;
  pMeshEnt e;
  apf::MeshIterator* ent_it = m->begin(dim);
  while ((e = m->iterate(ent_it)))
    if (((pGeomEnt)m->toModel(e))==g)  ents.push_back(e);
  m->end(ent_it);
}

// mesh entity functions
int pumi_ment_getdim(pMeshEnt e)
{
  return apf::getDimension(pumi::instance()->mesh, e);
}

int pumi_ment_getnumadj(pMeshEnt e, int target_dim)
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

void pumi_ment_getadj(pMeshEnt e, int target_dim, std::vector<pMeshEnt>& vecAdjEnt)
{
  int num_adj, ent_dim= apf::getDimension(pumi::instance()->mesh, e);
  if (ent_dim<target_dim) // upward
  {
    apf::Adjacent adjacent;
    pumi::instance()->mesh->getAdjacent(e,target_dim,adjacent);      
    std::copy(adjacent.begin(), adjacent.end(), vecAdjEnt.begin());
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

void pumi_ment_get2ndadj (pMeshEnt e, int bridge_dim, int target_dim, std::vector<pMeshEnt>& vecAdjEnt)
{
  apf::Adjacent adjacent;
  apf::getBridgeAdjacent(pumi::instance()->mesh, e, bridge_dim, target_dim, adjacent);
  std::copy(adjacent.begin(), adjacent.end(), vecAdjEnt.begin());
}


int pumi_ment_getlocalid(pMeshEnt e)
{
  return getMdsIndex(pumi::instance()->mesh, e);
}

pGeomEnt pumi_ment_getgeomclas(pMeshEnt e)
{
  return (pGeomEnt)(pumi::instance()->mesh->toModel(e));
}

pPartEnt pumi_ment_getptnclas(pMeshEnt e)
{  
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
  return NULL;
}

// owner part information
int pumi_ment_getownpid(pMeshEnt e)
{
  return pumi::instance()->mesh->getOwner(e);
}

pMeshEnt pumi_ment_getownent(pMeshEnt e)
{

  if (!(pumi::instance()->mesh->isShared(e))) // internal ent
    return e;

  int own_partid = pumi_ment_getownpid(e);
  if (own_partid==pumi_rank()) return e;

  return pumi_ment_getrmt(e, own_partid);
}

bool pumi_ment_isowned(pMeshEnt e)
{  
  if (pumi_ment_getownpid(e)==pumi_rank())
    return true;
  else
    return false;
}

// remote copy information
bool pumi_ment_isonbdry (pMeshEnt e)
{
  if (pumi::instance()->mesh->isShared(e)) // internal ent
    return true;
  else 
    return false;
}

int pumi_ment_getnumrmt (pMeshEnt e)
{
  if (!(pumi::instance()->mesh->isShared(e))) // internal ent
    return 0;
  apf::Copies remotes;
  pumi::instance()->mesh->getRemotes(e,remotes);
  return remotes.size();
}

void pumi_ment_getallrmt(pMeshEnt e, Copies& remotes)
{
  if (pumi::instance()->mesh->isShared(e))
    pumi::instance()->mesh->getRemotes(e,remotes);
}

pMeshEnt pumi_ment_getrmt(pMeshEnt& e, int pid)
{
  if (!(pumi::instance()->mesh->isShared(e))) // internal ent
    return NULL;

  Copies remotes;
  pumi::instance()->mesh->getRemotes(e,remotes);
  return remotes[pid];
}

void pumi_ment_setrmt(pMeshEnt e, int partID, pMeshEnt r)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

void pumi_ment_deletermt (pMeshEnt e, int partID)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

void pumi_ment_cleanrmt (pMeshEnt e)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

void pumi_ment_setptntopology (pMeshEnt e)
{
  if (!pumi_rank()) std::cout<<"[pumi error] "<<__func__<<" not supported\n";
}

int pumi_ment_getglobalid(pMeshEnt e)
{
  pTag tag = pumi::instance()->mesh->findTag("global_id");
  assert(tag);
  int global_id;
  pumi::instance()->mesh->getIntTag(e, tag, &global_id);
  return global_id;
  
}

void pumi_ment_getresidence(pMeshEnt e, std::vector<int>& resPartId)
{
  resPartId.clear();
  Copies remotes;
  pumi::instance()->mesh->getRemotes(e,remotes);
  APF_ITERATE(Copies,remotes,rit)
    resPartId.push_back(rit->first);
  resPartId.push_back(pumi_rank());
}

void pumi_ment_getclosureresidence(pMeshEnt e, std::vector<int>& resPartId)
{
  resPartId.clear();
  apf::Downward vertices;
  int nvtx = pumi::instance()->mesh->getDownward(e,0,vertices);
  pMeshEnt vtx;
  
  for (int i=0; i<nvtx; ++i)
  {
    vtx=vertices[i];
    Copies remotes;
    pumi::instance()->mesh->getRemotes(vtx,remotes);
    APF_ITERATE(Copies,remotes,rit)
      if (std::find(resPartId.begin(), resPartId.end(), rit->first)!=resPartId.end())
        resPartId.push_back(rit->first);
  }
  resPartId.push_back(pumi_rank());
}

// ghosting information
bool pumi_ment_isghost(pMeshEnt e)
{
  return pumi::instance()->mesh->isGhost(e);
}

bool pumi_ment_isghosted (pMeshEnt e)
{
  return pumi::instance()->mesh->isGhosted(e);
}

int pumi_ment_getnumghost (pMeshEnt e)
{
  if (!pumi::instance()->mesh->isGhosted(e)) return 0;
  apf::Copies ghosts;
  pumi::instance()->mesh->getGhosts(e,ghosts);
  return ghosts.size();
}

void pumi_ment_getallghost(pMeshEnt e, Copies& ghosts)
{
  if (pumi::instance()->mesh->isGhosted(e))
    pumi::instance()->mesh->getGhosts(e,ghosts);
}

pMeshEnt pumi_ment_getghost(pMeshEnt& e, int pid)
{
 if (!(pumi::instance()->mesh->isGhosted(e))) // internal ent
    return NULL;

  Copies ghosts;
  pumi::instance()->mesh->getGhosts(e,ghosts);
  return ghosts[pid];
}
