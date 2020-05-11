/****************************************************************************** 

  (c) 2004-2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include "gmi_mesh.h"
#include "gmi_null.h"
#include "gmi_analytic.h"
#include "pumi_iter.h"
#include "PCU.h"
#include <iostream>
#include <cstring>
#include <pcu_util.h>
#include <lionPrint.h>
#include "GenIterator.h"

gModel::gModel(gmi_model* model) : TagHolder() 
{
  g = model;
}

gModel::~gModel() {}

gEntity* gModel::getGeomEnt(int d, gmi_ent* ge)
{
  return allEntities.getGeomEnt(d, ge);
}

pGeom pumi_geom_load(const char* filename, const char* model_type, void (*geom_load_fp)(const char*))
{
  if (!strcmp(model_type,"null"))
  {
    gmi_register_null();
    return pumi_geom_load(gmi_load(".null"), model_type);
  }
  else if (!strcmp(model_type,"mesh"))
  {
    gmi_register_mesh();
    return pumi_geom_load(gmi_load(filename));
  }
  else if (!strcmp(model_type,"analytic")) 
    return pumi_geom_load(gmi_make_analytic(), model_type, filename, geom_load_fp);
  else
    if (!pumi_rank()) lion_eprint(1,"[PUMI ERROR] unsupported model type %s\n",model_type);
  
  return NULL;
}

pGeom pumi_geom_load(gmi_model* gm, const char* model_type, 
      const char* filename, void (*geom_load_fp)(const char*))
{
  double t0 = PCU_Time();
  if (!strcmp(model_type,"null"))
    pumi::instance()->model = new gModel(gm);
  else if (!strcmp(model_type,"mesh"))
  {
    pumi::instance()->model = new gModel(gm);
    pumi_geom_freeze(pumi::instance()->model);
  }
  else if (!strcmp(model_type,"analytic"))
  {
    pumi::instance()->model = new gModel(gm);
    if (geom_load_fp)
    {
      geom_load_fp(filename);
      pumi_geom_freeze(pumi::instance()->model);
    }
  }
  else
  {
    if (!pumi_rank()) lion_eprint(1,"[PUMI ERROR] unsupported model type %s\n",model_type);
    return NULL;
  }

  if (!PCU_Comm_Self() && filename)
    lion_oprint(1,"model %s loaded in %f seconds\n", filename, PCU_Time() - t0);

  return pumi::instance()->model;
}

void pumi_geom_delete(pGeom g)
{
  pTag id_tag=pumi_geom_findTag(g, "ID");

  for (int i=0; i<4; ++i) {
    std::vector<pGeomEnt> vge(g->size(i));
    for (pGeomIter gent_it = g->begin(i); gent_it!=g->end(i);++gent_it)
    {
      if (id_tag) pumi_gent_deleteTag(*gent_it, id_tag);
      vge.push_back(*gent_it);
    }
    for(size_t j=0; j<vge.size(); j++)
      delete vge[j];
  }
  pumi_geom_deleteTag(g, id_tag);
  delete g;
}

void pumi_geom_freeze(pGeom g)
{
  // loop over entities and fill the container
  for (int i=0; i<=3; ++i)
  {
    if (g->getGmi()->n[i]==g->size(i)) continue;
    PCU_ALWAYS_ASSERT(g->size(i)==0);
    gmi_iter* giter = gmi_begin(g->getGmi(), i);
    while(gmi_ent* gent = gmi_next(g->getGmi(), giter))
      g->add(i, new gEntity(gent));
    gmi_end(g->getGmi(), giter);
  }
}

void pumi_geom_createID(pGeom g) // ID starts from 1
{
  pTag id_tag = pumi_geom_createTag(g, "ID", PUMI_INT, 1);
  int id;
  for (int i=0; i<4; ++i)
  {
    id=1;
    for (pGeomIter gent_it = g->begin(i); gent_it!=g->end(i);++gent_it)
      pumi_gent_setIntTag(*gent_it, id_tag, id++);
  } 
}

pGeomEnt pumi_geom_findEnt(pGeom g, int d, int id)
{
  if (g->getGmi()->n[d]==0) 
    return NULL;
  gmi_ent* ge = gmi_find(g->getGmi(), d, id);
  if (ge) 
    return g->getGeomEnt(d, ge);
  else
    return NULL;
}

int pumi_geom_getNumEnt(pGeom g, int d)
{
  return g->getGmi()->n[d];
}

int pumi_giter_init (pGeom model, int type, gIter& iter)
{
  if (model->size(type)==0)
      return PUMI_FAILURE;

// FIXME: compilation error 
// /users/seol/develop/pumi/trunk/pumi/pumi_geom.cc:65:40: error: new initializer expression list treated as compound expression [-fpermissive] &processingNonFilter);
// /users/seol/develop/pumi/trunk/pumi/pumi_geom.cc:65:40: error: cannot convert 'void (*)(mPartEntityContainer::iter&, mPartEntityContainer::iter&, void*, int, int) {aka void (*)(ListIterator<gEntity>&, ListIterator<gEntity>&, void*, int, int)}' to 'gIter {aka GenIterator<ListIterator<gEntity>, gEntity>*}' in initialization

//int dim=model->size(3)?3:2;
//  iter = new gIter(model->begin(type), model->end(type), dim, type, (void*)model, 
//                   &processingNonFilter);
  if (iter->end())
    return PUMI_FAILURE;
  return PUMI_SUCCESS;
}

/* get next element of the mesh iterator */
int pumi_giter_getNext(gIter iter, pGeomEnt& ent)
{
   if(iter->end())
     return PUMI_FAILURE;
   else
     ent = **iter;
   iter->next();
   return PUMI_SUCCESS;
}

/* Check if the iterator has reached its end */
bool pumi_giter_isEnd(gIter iter)
{
  if(iter->end())
    return true;
  else
    return false;
}

void pumi_giter_delete(gIter iter)
{
  if (iter) delete iter;
  iter=NULL;
}

void pumi_giter_reset(gIter iter)
{
  iter->reset();
}

void pumi_geom_print (pGeom g, bool print_ent)
{
  if (PCU_Comm_Self()) return;
  std::cout<<"\n=== model entity and tag info === \n";
  std::cout<<"# global geom ent: v "<<g->size(0)<<", e "
           <<g->size(1)<<", f "<<g->size(2)<<", r "<<g->size(3)<<"\n";

  std::vector<pTag> tags;
  pumi_geom_getTag (g, tags);
  int n = (int) tags.size();
  for (int i = 0; i < n; ++i) 
    std::cout<<"tag "<<i<<": \""<< Tag_GetName(tags[i])<<"\", type "
             << pumi_tag_getType(tags[i])<<", size "<< pumi_tag_getSize(tags[i])<<"\n";

  // print model entities
  if (!print_ent) return;

  const char* name[4] = {"vtx", "edge", "face", "rgn"};
  for (int d=0; d<4; ++d)
  {
    for (pGeomIter gent_it = g->begin(d); gent_it!=g->end(d);++gent_it)
    {
      std::cout<<"geom "<<name[d]<<" "<<pumi_gent_getID(*gent_it)<<": ";
      if (d) 
      {
        for (int adj_d=0; adj_d<d; ++adj_d)
        {
          std::vector<pGeomEnt> adj;
          pumi_gent_getAdj (*gent_it, adj_d, adj);
          if (adj.size()) 
          {
            if (adj_d) std::cout<<", ";
            std::cout<<name[adj_d];
            for (std::vector<pGeomEnt>::iterator git=adj.begin(); git!=adj.end(); ++git)
              std::cout<<" "<<pumi_gent_getID(*git);
          }        
        }
      }
      std::cout<<"\n";
    }
  }
  std::cout<<"\n";
}
