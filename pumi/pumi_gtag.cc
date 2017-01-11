/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <assert.h>
#include "GenTag.h"
#include <iostream>

//**********************************************************
int pumi_tag_getType (const pTag tag)
//**********************************************************
{
  return Tag_GetType(tag);
}

//**********************************************************
void pumi_tag_getName (const pTag tag, const char** name)
//**********************************************************
{
  *name = Tag_GetName(tag);
}

//**********************************************************
int pumi_tag_getSize (const pTag tag)
//**********************************************************
{
  return Tag_GetSize(tag);
}

//**********************************************************
int pumi_tag_getByte (const pTag tag)
//**********************************************************
{
  return tag->getBytes();
}

//**********************************************************
pTag pumi_geom_createTag (pGeom m, const char* tag_name, int tag_type, int tag_size)
//**********************************************************
{
  pTag t = pumi_geom_findTag (m, tag_name);
  if (t)  return t;

  return TagHolder_CreateTag(static_cast<pTagHolder>(m), tag_name, tag_type, tag_size); 
}

//**********************************************************
void pumi_geom_deleteTag (pGeom g, pTag tag, bool force_delete)
//**********************************************************
{ 
  if (!pumi_geom_hasTag (g, tag)) return;

  if (force_delete)
  {
    for (int i=0; i<4; ++i)
      for (pGeomIter git = g->begin(i); git!=g->end(i);++git)
        pumi_gent_deleteTag (*git, tag);
  }

  // if force_delete==0 and tag handle is still in use, pumi_geom_deleteTag crashes
  TagHolder_DelTag(static_cast<pTagHolder>(g), tag); 
}

//**********************************************************
pTag pumi_geom_findTag (pGeom m, const char* name)
//**********************************************************
{
  pTag tag;
  if (!TagHolder_FindTag (static_cast<pTagHolder>(m), name, tag))
    return tag;
  else
    return NULL;
}

//**********************************************************
bool pumi_geom_hasTag (pGeom m, const pTag tag)
//**********************************************************
{ 
  int exist;
  assert(!TagHolder_HasTag (static_cast<pTagHolder>(m), tag, &exist));
  if (exist) 
    return true;
  else
    return false;
}

//**********************************************************
void pumi_geom_getTag (pGeom m, std::vector<pTag>& tags)
//**********************************************************
{
  TagHolder_GetTag(static_cast<pTagHolder>(m), tags);
}

// entity tagging
//**********************************************************
bool pumi_gent_hasTag (pGeomEnt ent, pTag tag)
//**********************************************************
{
  return Taggable_HasTag(static_cast<pTaggable>(ent), tag);
}

//**********************************************************
void pumi_gent_deleteTag (pGeomEnt ent, pTag tag)
//**********************************************************
{
  Taggable_DelTag(static_cast<pTaggable>(ent), tag);
}

//**********************************************************
void pumi_gent_getTag (pGeomEnt ent, std::vector<pTag>& tags)
//**********************************************************
{
  Taggable_GetTag (static_cast<pTaggable>(ent), tags);
}

//**********************************************************
void pumi_gent_setStringTag(pGeomEnt ent, pTag tag, const char* s)
//**********************************************************
{
  ent->setTagString(tag,s);
}

//**********************************************************
void pumi_gent_getStringTag(pGeomEnt ent, pTag tag, const char*& s)
//**********************************************************
{
  s = ent->getTagString(tag);
}

//**********************************************************
void pumi_gent_setPtrTag (pGeomEnt ent, pTag tag, void* data)
//**********************************************************
{ 
  assert(!Taggable_SetData<void*>(static_cast<pTaggable>(ent), tag, &data));
 }

//**********************************************************
void pumi_gent_getPtrTag (pGeomEnt ent, pTag tag, void** data)
//**********************************************************
{ 
  assert(!Taggable_GetData<void*>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_setIntTag (pGeomEnt ent, pTag tag, const int data)
//**********************************************************
{
  assert(!Taggable_SetData<int>(static_cast<pTaggable>(ent), tag, &data));
}

//**********************************************************
void pumi_gent_getIntTag (pGeomEnt ent, pTag tag, int* data)
//**********************************************************
{
  assert(!Taggable_GetData<int>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_setLongTag (pGeomEnt ent, pTag tag, const long data)
//**********************************************************
{
  assert(!Taggable_SetData<long>(static_cast<pTaggable>(ent), tag, &data));
}

//**********************************************************
void pumi_gent_getLongTag (pGeomEnt ent, pTag tag, long* data)
//**********************************************************
{
  assert(!Taggable_GetData<long>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_setDblTag (pGeomEnt ent, pTag tag, const double data)
//**********************************************************
{
  assert(!Taggable_SetData<double>(static_cast<pTaggable>(ent), tag, &data));
}

//**********************************************************
void pumi_gent_getDblTag (pGeomEnt ent, pTag tag, double* data)
//**********************************************************
{
  assert(!Taggable_GetData<double>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_setEntTag (pGeomEnt ent, pTag tag, const pGeomEnt data)
//**********************************************************
{
  assert(!Taggable_SetData<pGeomEnt>(static_cast<pTaggable>(ent), tag, &data));
}

//**********************************************************
void pumi_gent_getEntTag (pGeomEnt ent, pTag tag, pGeomEnt *data)
//**********************************************************
{
  assert(!Taggable_GetData<pGeomEnt>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_setPtrArrTag (pGeomEnt ent, pTag tag, void* const* data)
//**********************************************************
{
  assert(!Taggable_SetData<void*>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_getPtrArrTag (pGeomEnt ent, pTag tag, void** data)
//**********************************************************
{
  assert(!Taggable_GetData<void*>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_setIntArrTag (pGeomEnt ent, pTag tag, const int* data)
//**********************************************************
{
  assert(!Taggable_SetData<int>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_getIntArrTag (pGeomEnt ent, pTag tag, int** data, int* data_size)
//**********************************************************
{
  *data_size=Tag_GetSize(tag);
  assert(!Taggable_GetData<int>(static_cast<pTaggable>(ent), tag, *data));
}

//**********************************************************
void pumi_gent_setLongArrTag (pGeomEnt ent, pTag tag, const long* data)
//**********************************************************
{
  assert(!Taggable_SetData<long>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_getLongArrTag (pGeomEnt ent, pTag tag, long** data, int* data_size)
//**********************************************************
{
  *data_size=Tag_GetSize(tag);
  assert(!Taggable_GetData<long>(static_cast<pTaggable>(ent), tag, *data));
}

//**********************************************************
void pumi_gent_setDblArrTag (pGeomEnt ent, pTag tag, const double* data)
//**********************************************************
{
  assert(!Taggable_SetData<double>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_getDblArrTag (pGeomEnt ent, pTag tag, double** data, int* data_size)
//**********************************************************
{
  *data_size=Tag_GetSize(tag);
  assert(!Taggable_GetData<double>(static_cast<pTaggable>(ent), tag, *data));
}

//**********************************************************
void pumi_gent_setEntArrTag (pGeomEnt ent, pTag tag, const pGeomEnt* data)
//**********************************************************
{
  assert(!Taggable_SetData<pGeomEnt>(static_cast<pTaggable>(ent), tag, data));
}

//**********************************************************
void pumi_gent_getEntArrTag (pGeomEnt ent, pTag tag, pGeomEnt** data, int* data_size)
//**********************************************************
{
  *data_size=Tag_GetSize(tag);
  assert(!Taggable_GetData<pGeomEnt>(static_cast<pTaggable>(ent), tag, *data));
}
