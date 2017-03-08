/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfTagData.h"
#include "apfShape.h"

#include <pcu_util.h>

namespace apf {

void TagData::init(
    const char* name,
    Mesh* m,
    FieldShape* s,
    TagMaker* mk,
    int c)
{
  mesh = m;
  shape = s;
  maker = mk;
  createTags(name,c);
}

TagData::~TagData()
{
  detachTags();
  destroyTags();
}

bool TagData::hasEntity(MeshEntity* e)
{
  MeshTag* tag = getTag(e);
  if (!tag) return false;
  return mesh->hasTag(e,tag);
}

void TagData::removeEntity(MeshEntity* e)
{
  MeshTag* tag = this->getTag(e);
  if (!tag) return;
  mesh->removeTag(e,this->getTag(e));
}

MeshTag* TagData::getTag(MeshEntity* e)
{
  return tags[mesh->getType(e)];
}

MeshTag* TagData::makeOrFindTag(const char* name, int size)
{
  MeshTag* tag = mesh->findTag(name);
  if (tag) return tag;
  return maker->make(mesh,name,size);
}

static const char* typePostfix[Mesh::TYPES] =
{"ver","edg","tri","qua","tet","hex","pri","pyr"};

void TagData::createTags(const char* name, int components)
{
  PCU_ALWAYS_ASSERT(name);
  PCU_ALWAYS_ASSERT(shape);
  for (int type=Mesh::VERTEX; type < Mesh::TYPES; ++type)
  {
    int n = shape->countNodesOn(type);
    if (n)
    {
      std::string tagName(name);
      tagName += '_';
      tagName += typePostfix[type];
      tags[type] = makeOrFindTag(tagName.c_str(),n*components);
    }
    else
      tags[type] = 0;
  }
}

void TagData::detachTags()
{
  for (int dim=0; dim <= 3; ++dim)
    if (shape->hasNodesIn(dim))
    {
      MeshEntity* entity;
      MeshIterator* entities = mesh->begin(dim);
      while ((entity = mesh->iterate(entities)))
        this->removeEntity(entity);
      mesh->end(entities);
    }
}

void TagData::destroyTags()
{
  for (int type=Mesh::VERTEX; type < Mesh::TYPES; ++type)
    if (tags[type])
      mesh->destroyTag(tags[type]);
}

void TagData::rename(const char* newName)
{
  for (int type=Mesh::VERTEX; type < Mesh::TYPES; ++type)
    if (tags[type]) {
      std::string newTagName(newName);
      newTagName += '_';
      newTagName += typePostfix[type];
      mesh->renameTag(tags[type], newTagName.c_str());
    }
}

}
