/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include "gmi_mesh.h"
#include "gmi_null.h"
#include "PCU.h"
#include <iostream>
#include <cstring>
pGeom pumi_geom_create(const char* filename, const char* model_type)
{
  if (!strcmp(model_type,"mesh"))
  {
    gmi_register_mesh();
    pumi::instance()->model = gmi_load(filename);
  }

  else if (!strcmp(model_type,"null"))
  {
    gmi_register_null();
    pumi::instance()->model = gmi_load(".null");
  }
  else
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid model type "<<model_type<<"\n";

  return pumi::instance()->model;
}
