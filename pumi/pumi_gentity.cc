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
#include <vector>
#include <algorithm>
#include <iostream>
#include <pcu_util.h>

using std::vector;

gEntity::gEntity(gmi_ent* ent) : Taggable()
{ 
  e = ent; 
}

gEntity::~gEntity()
{}
 
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

void pumi_gent_getAdj (pGeomEnt ge, int target_dim, std::vector<pGeomEnt>& adj_ents_vec)
{
  pGeom g = pumi::instance()->model;
  int ent_dim = gmi_dim(g->getGmi(), ge->getGmi());
  if (ent_dim==target_dim) return;

  std::set<pGeomEnt> adj_ents;
  std::set<pGeomEnt> ents;
  ents.insert(ge);

  if (abs(ent_dim-target_dim)==1) 
    get_one_level_adj (g, ents, target_dim, adj_ents);
  else
  {
    switch(ent_dim)
    {
      case 0: if (target_dim==2) 
            {
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, ents, 1, edges);
              get_one_level_adj (g, edges, 2, adj_ents);
            }  
            if (target_dim==3)
            {
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, ents, 1, edges);
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, edges, 2, faces);
              get_one_level_adj (g, faces, 3, adj_ents);
            }  
            break;
      case 1: if (target_dim==3)
            {
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, ents, 2, faces);
              get_one_level_adj (g, faces, 3, adj_ents);
            }
            break;
      case 2: if (target_dim==0)
            {
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, ents, 1, edges);
              get_one_level_adj (g, edges, 0, adj_ents);
            }  
            break;
      case 3: if (target_dim==0)
            {
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, ents, 2, faces);
              std::set<pGeomEnt> edges;
              get_one_level_adj (g, faces, 1, edges);
              get_one_level_adj (g, edges, 0, adj_ents);
            }  
            if (target_dim==1)
            {
              std::set<pGeomEnt> faces;
              get_one_level_adj (g, ents, 2, faces);
              get_one_level_adj (g, faces, 1, adj_ents);
            }  
            break;
      default: break;
    }
  }
  for (std::set<pGeomEnt>::iterator it=adj_ents.begin(); it!=adj_ents.end(); ++it)
    adj_ents_vec.push_back(*it);
}

int pumi_gent_getNumAdj (pGeomEnt ge, int target_dim)
{
  std::vector<pGeomEnt> adj_ents;
  pumi_gent_getAdj (ge, target_dim, adj_ents);
  return (int)(adj_ents.size());
}
