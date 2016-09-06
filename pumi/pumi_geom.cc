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
#include "PCU.h"
#include <iostream>
#include <cstring>
#include <assert.h>

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
