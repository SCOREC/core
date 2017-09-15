/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <pcu_util.h>

using std::vector;

gEntity::gEntity(gmi_ent* ent) : Taggable()
{ 
  e = ent; 
}

gEntity::~gEntity() {}
 
int pumi_gent_getDim(pGeomEnt ge)
{
  return gmi_dim(pumi::instance()->model->getGmi(), ge->getGmi());
}

int pumi_gent_getID(pGeomEnt ge)
{
  return gmi_tag(pumi::instance()->model->getGmi(), ge->getGmi());
}

void pumi_gent_getRevClas (pGeomEnt ge, std::vector<pMeshEnt>& ents)
{
  PCU_ALWAYS_ASSERT(!ents.size());
  int dim=gmi_dim(pumi::instance()->model->getGmi(), ge->getGmi());
  pMesh m = pumi::instance()->mesh;
  pMeshEnt e;
  apf::MeshIterator* ent_it = m->begin(dim);
  while ((e = m->iterate(ent_it)))
    if (((gmi_ent*)m->toModel(e))==ge->getGmi())  ents.push_back(e);
  m->end(ent_it);
}

void get_one_level_adj (pGeom g, std::set<pGeomEnt>& ents, 
                        int dim, std::set<pGeomEnt>& adj_ents)
{
  for (std::set<pGeomEnt>::iterator it=ents.begin(); it!=ents.end(); ++it)
  {
    gmi_set* g_adj = gmi_adjacent(g->getGmi(), (*it)->getGmi(), dim);
    for (int i=0; i<g_adj->n; ++i)
      adj_ents.insert(g->getGeomEnt(dim, g_adj->e[i]));
  }
}

void getAdjacent (pGeom g, pGeomEnt ge, int target_dim, std::set<pGeomEnt>& result)
{
  int ent_dim = gmi_dim(g->getGmi(), ge->getGmi());
  PCU_ALWAYS_ASSERT(ent_dim != target_dim);
  std::set<pGeomEnt> ents;
  ents.insert(ge);

  if (abs(ent_dim-target_dim)==1) 
  {
    get_one_level_adj (g, ents, target_dim, result);
    return;
  }

  //  if not (abs(ent_dim-target_dim)==1) 
  switch(ent_dim)
  {
    case 0: if (target_dim==2) 
            {
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, ents, 1, edges);
              get_one_level_adj (g, edges, 2, result);
            }  
            if (target_dim==3)
            {
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, ents, 1, edges);
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, edges, 2, faces);
              get_one_level_adj (g, faces, 3, result);
            }  
            break;
    case 1: if (target_dim==3)
            {
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, ents, 2, faces);
              get_one_level_adj (g, faces, 3, result);
            }
            break;
    case 2: if (target_dim==0)
            {
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, ents, 1, edges);
              get_one_level_adj (g, edges, 0, result);
            }  
            break;
    case 3: if (target_dim==0)
            {
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, ents, 2, faces);
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, faces, 1, edges);
              get_one_level_adj (g, edges, 0, result);
            }  
            if (target_dim==1)
            {
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, ents, 2, faces);
              get_one_level_adj (g, faces, 1, result);
            }  
            break;
    default: break;
  } // switch

}

void pumi_gent_getAdj (pGeomEnt ge, int target_dim, std::vector<pGeomEnt>& result_vec)
{
  std::set<pGeomEnt> result;
  getAdjacent (pumi::instance()->model, ge, target_dim, result);
 
  for (std::set<pGeomEnt>::iterator git=result.begin(); git!=result.end(); ++git)
    result_vec.push_back(*git);
}

int pumi_gent_getNumAdj (pGeomEnt ge, int target_dim)
{
  std::set<pGeomEnt> result;
  getAdjacent (pumi::instance()->model, ge, target_dim, result);
  return (int)(result.size());
}

void pumi_gent_get2ndAdj (pGeomEnt ge, 
     int bridgeDimension, int targetDimension, std::vector<pGeomEnt>& result_vec)
{
  pGeom g = pumi::instance()->model;
  PCU_ALWAYS_ASSERT(targetDimension != bridgeDimension);
  std::set<pGeomEnt> result;
  std::set<pGeomEnt> bridges;
  getAdjacent(g, ge, bridgeDimension, bridges);
  for (std::set<pGeomEnt>::iterator bit=bridges.begin(); bit!=bridges.end(); ++bit)
  {
    std::set<pGeomEnt> targets;
    getAdjacent(g, *bit, targetDimension, targets);
    for (std::set<pGeomEnt>::iterator tit=targets.begin(); tit!=targets.end(); ++tit)
      result.insert(*tit);
  }
  result.erase(ge);
  for (std::set<pGeomEnt>::iterator git=result.begin(); git!=result.end(); ++git)
    result_vec.push_back(*git);
}
