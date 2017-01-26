/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
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
#include <assert.h>
#include "GenIterator.h"

gModel::gModel(gmi_model* model) : TagHolder() 
{
  g = model;
}

gModel::~gModel()
{
  delete g;
}

gEntity* gModel::getGeomEnt(gmi_ent* ge)
{
  return allEntities.getGeomEnt(gmi_dim(g, ge), ge);
}

pGeom pumi_geom_load(const char* filename, const char* model_type, void (*geom_load_fp)(const char*))
{
  if (!strcmp(model_type,"null"))
  {
    gmi_register_null();
    pumi::instance()->model = new gModel(gmi_load(".null"));
  }
  else if (!strcmp(model_type,"mesh"))
  {
    gmi_register_mesh();
    pumi::instance()->model = new gModel(gmi_load(filename));
    pumi_geom_freeze(pumi::instance()->model);
  }
  else if (!strcmp(model_type,"analytic")) 
  {
    pumi::instance()->model = new gModel(gmi_make_analytic());
    if (geom_load_fp)
    {
      geom_load_fp(filename);
      pumi_geom_freeze(pumi::instance()->model);
    }
  }
  else
  {
    if (!pumi_rank()) std::cerr<<"[PUMI ERROR] unsupported model type "<<model_type<<"\n";
    return NULL;
  }

  return pumi::instance()->model;
}

void pumi_geom_freeze(pGeom g)
{
  // loop over entities and fill the container
  for (int i=0; i<=3; ++i)
  {
    if (g->getGmi()->n[i]==g->size(i)) continue;
    assert(g->size(i)==0);
    gmi_iter* giter = gmi_begin(g->getGmi(), i);
    while(gmi_ent* gent = gmi_next(g->getGmi(), giter))
      g->add(i, new gEntity(gent));
    gmi_end(g->getGmi(), giter);
  }
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
