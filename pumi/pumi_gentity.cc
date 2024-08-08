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
  pTag id_tag=pumi_geom_findTag(pumi::instance()->model, "ID");
  if (id_tag) 
  { 
    int int_data;
    pumi_gent_getIntTag (ge, id_tag, &int_data);
    return int_data;
  }
  else
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

void get_one_level_adj (gmi_model* g, std::set<gmi_ent*>& ents, 
                        int dim, std::set<gmi_ent*>& adj_ents)
{
  for (std::set<gmi_ent*>::iterator it=ents.begin(); it!=ents.end(); ++it)
  {
    gmi_set* g_adj = gmi_adjacent(g, *it, dim);
    for (int i=0; i<g_adj->n; ++i)
      adj_ents.insert(g_adj->e[i]);
    gmi_free_set(g_adj);
  }
}

void gmi_getAdj(gmi_model* g, gmi_ent* ge, int tgt_dim, std::set<gmi_ent*>& result)
{
  int ent_dim = gmi_dim(g, ge);
  PCU_ALWAYS_ASSERT(ent_dim != tgt_dim);
  std::set<gmi_ent*> ents;
  ents.insert(ge);

  if (abs(ent_dim-tgt_dim)==1) 
  {
    get_one_level_adj (g, ents, tgt_dim, result);
    return;
  }

  //  if not (abs(ent_dim-tgt_dim)==1) 
  switch(ent_dim)
  {
    case 0: if (tgt_dim==2) 
            {
              std::set<gmi_ent*> edges;
              get_one_level_adj (g, ents, 1, edges);
              get_one_level_adj (g, edges, 2, result);
            }  
            if (tgt_dim==3)
            {
              std::set<gmi_ent*> edges;
              get_one_level_adj (g, ents, 1, edges);
              std::set<gmi_ent*> faces;
              get_one_level_adj (g, edges, 2, faces);
              get_one_level_adj (g, faces, 3, result);
            }  
            break;
    case 1: if (tgt_dim==3)
            {
              std::set<gmi_ent*> faces;
              get_one_level_adj (g, ents, 2, faces);
              get_one_level_adj (g, faces, 3, result);
            }
            break;
    case 2: if (tgt_dim==0)
            {
              std::set<gmi_ent*> edges;
              get_one_level_adj (g, ents, 1, edges);
              get_one_level_adj (g, edges, 0, result);
            }  
            break;
    case 3: if (tgt_dim==0)
            {
              std::set<gmi_ent*> faces;
              get_one_level_adj (g, ents, 2, faces);
              std::set<gmi_ent*> edges;
              get_one_level_adj (g, faces, 1, edges);
              get_one_level_adj (g, edges, 0, result);
            }  
            if (tgt_dim==1)
            {
              std::set<gmi_ent*> faces;
              get_one_level_adj (g, ents, 2, faces);
              get_one_level_adj (g, faces, 1, result);
            }  
            break;
    default: break;
  } // switch
}

void gmi_get2ndAdj (gmi_model* g, gmi_ent* ge, 
     int brg_dim, int tgt_dim, std::set<gmi_ent*>& result)
{
  PCU_ALWAYS_ASSERT(tgt_dim != brg_dim && result.empty());
  std::set<gmi_ent*> bridges;
  gmi_getAdj(g, ge, brg_dim, bridges);
  for (std::set<gmi_ent*>::iterator bit=bridges.begin(); bit!=bridges.end(); ++bit)
  {
    std::set<gmi_ent*> targets;
    gmi_getAdj(g, *bit, tgt_dim, targets);
    for (std::set<gmi_ent*>::iterator tit=targets.begin(); tit!=targets.end(); ++tit)
      result.insert(*tit);
  }
  result.erase(ge);
}


void pumi_gent_getAdj (pGeomEnt ge, int tgt_dim, std::vector<pGeomEnt>& result_vec)
{
  std::set<gmi_ent*> result;
  gmi_getAdj(pumi::instance()->model->getGmi(), ge->getGmi(), tgt_dim, result);
 
  for (std::set<gmi_ent*>::iterator git=result.begin(); git!=result.end(); ++git)
    result_vec.push_back(pumi::instance()->model->getGeomEnt(tgt_dim, *git));
}

int pumi_gent_getNumAdj (pGeomEnt ge, int tgt_dim)
{
  std::set<gmi_ent*> result;
  gmi_getAdj(pumi::instance()->model->getGmi(), ge->getGmi(), tgt_dim, result);
  return (int)(result.size());
}

void pumi_gent_get2ndAdj (pGeomEnt ge, 
     int brg_dim, int tgt_dim, std::vector<pGeomEnt>& result_vec)
{
  std::set<gmi_ent*> result;
  gmi_get2ndAdj (pumi::instance()->model->getGmi(), ge->getGmi(), brg_dim, tgt_dim, result);

  for (std::set<gmi_ent*>::iterator git=result.begin(); git!=result.end(); ++git)
    result_vec.push_back(pumi::instance()->model->getGeomEnt(tgt_dim, *git));
}
