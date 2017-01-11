/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"

pMeshTag pumi_mesh_createIntTag(pMesh m, const char* name, int size)
{ 
  return m->createIntTag(name, size);
}

pMeshTag pumi_mesh_createLongTag(pMesh m, const char* name, int size)
{ 
  return m->createLongTag(name, size);
}

pMeshTag pumi_mesh_createDblTag(pMesh m, const char* name, int size)
{ 
  return m->createDoubleTag(name, size);
}

//*******************************************************
void pumi_mesh_deleteTag(pMesh m, pMeshTag t, bool force_delete)
//*******************************************************
{
  if (force_delete)
    for (int i=0; i<4; ++i)
      apf::removeTagFromDimension(m, t, i);

  m->destroyTag(t);
}

pMeshTag pumi_mesh_findTag(pMesh m, const char* name)
{
  return m->findTag(name); 
}

bool pumi_mesh_hasTag (pMesh m, const pMeshTag t)
{
  if (m->findTag(m->getTagName(t))) 
    return true;
  else
    return false;
}

void pumi_mesh_getTag(pMesh m, std::vector<pMeshTag> tag_vec)
{
  apf::DynamicArray<pMeshTag> tags;
  m->getTags(tags);
  for (size_t n = 0; n<tags.getSize();++n)
    tag_vec.push_back(tags[n]);
}

// tag management over mesh entity
void pumi_ment_deleteTag (pMeshEnt e, pMeshTag t)
{
  pumi::instance()->mesh->removeTag(e, t);
}

bool pumi_ment_hasTag (pMeshEnt e, pMeshTag t)
{
  return pumi::instance()->mesh->hasTag(e, t);
}

// set/get tag data to/from mesh entity
void pumi_ment_setIntTag(pMeshEnt e, pMeshTag t, int const* data)
{
  pumi::instance()->mesh->setIntTag(e, t, data);
}

void pumi_ment_getIntTag(pMeshEnt e, pMeshTag t, int* data)
{
  pumi::instance()->mesh->getIntTag(e, t, data);
}


void pumi_ment_setLongTag(pMeshEnt e, pMeshTag t, long const* data)
{
  pumi::instance()->mesh->setLongTag(e, t, data);
}

void pumi_ment_getLongTag(pMeshEnt e, pMeshTag t, long* data)
{
  pumi::instance()->mesh->getLongTag(e, t, data);
}

void pumi_ment_setDblTag(pMeshEnt e, pMeshTag t, double const* data)
{
  pumi::instance()->mesh->setDoubleTag(e, t, data);
}

void pumi_ment_getDblTag(pMeshEnt e, pMeshTag t, double* data)
{
  pumi::instance()->mesh->getDoubleTag(e, t, data);
}

