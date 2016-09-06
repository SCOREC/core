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
#include <assert.h>

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
  assert(!ents.size());
  int dim=gmi_dim(pumi::instance()->model->getGmi(), ge->getGmi());
  pMesh m = pumi::instance()->mesh;
  pMeshEnt e;
  apf::MeshIterator* ent_it = m->begin(dim);
  while ((e = m->iterate(ent_it)))
    if (((gmi_ent*)m->toModel(e))==ge->getGmi())  ents.push_back(e);
  m->end(ent_it);
}
