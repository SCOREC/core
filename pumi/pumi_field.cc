/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include "apf.h"

pField pumi_field_create(pMesh m, const char* name, int num_dof_per_elm)
{
  return apf::createPackedField(m, name, num_dof_per_elm);
}

int pumi_field_getSize(pField f)
{  
  return apf::countComponents(f); 
}

int pumi_field_getType(pField f)
{ 
  return apf::getValueType(f); 
}

std::string pumi_field_getName(pField f)
{ 
  return apf::getName(f); 
}

void pumi_field_delete(pField f)
{
  apf::destroyField(f);
}

void pumi_field_synchronize(apf::Field* f)
{  
  apf::synchronize(f, getSharing(getMesh(f)));
}

void pumi_field_accumulate(apf::Field* f)
{  
  apf::accumulate(f, getSharing(getMesh(f)));
}

void pumi_field_freeze(apf::Field* f)
{  
  apf::freeze(f); 
}

void pumi_field_unfreeze(apf::Field* f)
{  
  apf::unfreeze(f); 
}

//*******************************************************
void pumi_ment_getField (pMeshEnt e, pField f, double* dof_data)
//*******************************************************
{
  apf::getComponents(f, e, 0, dof_data);
}    

void pumi_ment_setField (pMeshEnt e, pField f, double* dof_data)
{
  apf::setComponents(f, e, 0, dof_data);
}
