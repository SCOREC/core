/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "mPartEntityContainer.h"
#include <pcu_util.h>

mPartEntityContainer::mPartEntityContainer()
{
}

mPartEntityContainer::~mPartEntityContainer()
{
}

mPartEntityContainer::iter mPartEntityContainer::begin(int what)
{
  PCU_ALWAYS_ASSERT(what >= 0);
  PCU_ALWAYS_ASSERT(what < _DIMS_);
  return gEntities[what].begin<gEntity>();
}

mPartEntityContainer::iter mPartEntityContainer::end(int what)
{
  PCU_ALWAYS_ASSERT(what >= 0);
  PCU_ALWAYS_ASSERT(what < _DIMS_);
  return gEntities[what].end<gEntity>();
}

void mPartEntityContainer::add(int d, gEntity* e)
{ 
  gEntities[d].push_back(e);
  gEntities_map[d].insert(std::map<gmi_ent*, gEntity*>::value_type(e->getGmi(), e));
}

void mPartEntityContainer::del(int d, gEntity* e)
{ 
  gEntities[d].remove(e);
  gEntities_map[d].erase(std::map<gmi_ent*, gEntity*>::key_type(e->getGmi()));
}

gEntity* mPartEntityContainer::getGeomEnt(int d, gmi_ent* e)
{
  return gEntities_map[d][e];
}


int mPartEntityContainer::size(int what) const 
{
  PCU_ALWAYS_ASSERT(what >= 0);
  PCU_ALWAYS_ASSERT(what < _DIMS_);
  return gEntities[what].size();
} 

