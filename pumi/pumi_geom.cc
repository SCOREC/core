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
pGeom pumi_geom_load(const char* filename, const char* model_type)
{
  assert(!strcmp(model_type,"null") || !strcmp(model_type,"mesh") || !strcmp(model_type,"analytic"));
  if (!strcmp(model_type,"null"))
  {
    gmi_register_null();
    pumi::instance()->model = new gModel(gmi_load(".null"));
  }
  else
  {
    if (!strcmp(model_type,"mesh"))
      gmi_register_mesh();
    pumi::instance()->model = new gModel(gmi_load(filename));
  }
  
  // loop over entities and fill the container
  for (int i=0; i<=3; ++i)
  {
    gmi_iter* giter = gmi_begin(pumi::instance()->model->getGmi(), i);
    while(gmi_ent* gent = gmi_next(pumi::instance()->model->getGmi(), giter))
      pumi::instance()->model->add(i, new gEntity(gent));
    gmi_end(pumi::instance()->model->getGmi(), giter);
  }
  return pumi::instance()->model;
}

inline void processingNonFilter(mPartEntityContainer::iter& it_begin, mPartEntityContainer::iter& it_end, void* ptr, int type, int topo)
{}

int pumi_giter_init (pGeom model, int type, gIter& iter)
{
  if (model->size(type)==0)
      return PUMI_FAILURE;
  int dim=model->size(3)?3:2;

// FIXME: compilation error 
// /users/seol/develop/pumi/trunk/pumi/pumi_geom.cc:65:40: error: new initializer expression list treated as compound expression [-fpermissive] &processingNonFilter);
// /users/seol/develop/pumi/trunk/pumi/pumi_geom.cc:65:40: error: cannot convert 'void (*)(mPartEntityContainer::iter&, mPartEntityContainer::iter&, void*, int, int) {aka void (*)(ListIterator<gEntity>&, ListIterator<gEntity>&, void*, int, int)}' to 'gIter {aka GenIterator<ListIterator<gEntity>, gEntity>*}' in initialization

//  iter = new gIter(model->beginall(type), model->endall(type), dim, type, (void*)model, 
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
