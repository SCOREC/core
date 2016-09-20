/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "GenTag.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

size_t const TagHandle::typeSizes[PUMI_TAGTYPES] = 
{sizeof(char)    //PUMI_BYTE
,sizeof(int)     //PUMI_INT
,sizeof(double)  //PUMI_DBL
,sizeof(void*)   //PUMI_ENT
,sizeof(void*)   //PUMI_SET
,sizeof(void*)   //PUMI_PTR
,0               //PUMI_STR
,sizeof(long)    //PUMI_LONG
};

TagHandle::TagHandle(
    std::string const& in_name,
    const int in_type,
    const int in_size)
{
  tag_name = in_name;
  tag_type = in_type;
  tag_size = in_size;
  bytes = typeSizes[in_type]*in_size;
}

TagHandle::~TagHandle()
{
}

bool TagHandle::operator<(const TagHandle& other) const
{
  return tag_name < other.tag_name;
}

Taggable::Taggable()
{
  container = 0;
  size = 0;
}

Taggable::~Taggable()
{
  clearTagData();
}

typedef Taggable::Entry Entry;

static void removeEntry(Entry*& container, int e, int& size)
{
  Entry* container2 = new Entry[size-1];
  int i;
  int j=0;
  for (i=0; i < e; ++i)
    container2[j++] = container[i];
  for (i=e+1; i < size; ++i)
    container2[j++] = container[i];
  delete [] container;
  container = container2;
  --size;
}

static Entry* addEntry(
    Entry*& container,
    TagHandle* handle,
    int& size)
{
  Entry* container2 = new Entry[size+1];
  for (int i=0; i < size; ++i)
    container2[i] = container[i];
  container2[size].handle = handle;
  delete [] container;
  container = container2;
  Entry* e = container + size;
  ++size;
  return e;
}

static int findEntry(Entry* container, TagHandle* handle, int size)
{
  for (int i=0; i < size; ++i)
    if (container[i].handle == handle)
      return i;
  return -1;
}

static void freeEntry(Entry* e)
{
  TagHandle* handle = e->handle;
  if ((handle->getBytes() > sizeof(void*))
    ||(handle->getType() == PUMI_STR))
    free(e->data);
}

void Taggable::deleteTagData(TagHandle* tag)
{
  int const i = findEntry(container,tag,size);
  if (i==-1)
    return;
  freeEntry(container+i);
  removeEntry(container,i,size);
}

void Taggable::clearTagData()
{
  for (int i=0; i < size; ++i)
    freeEntry(container+i);
  delete [] container;
}

bool Taggable::hasTagData(TagHandle* tag)
{
  return findEntry(container,tag,size) != -1;
}

bool Taggable::getTagData(TagHandle* tag, void* data)
{
  int const i = findEntry(container,tag,size);
  if (i==-1)
    return false;
  Entry* e = container + i;
  size_t const bytes = e->handle->getBytes();
  if (bytes > sizeof(void*))
    memcpy(data,e->data,bytes);
  else
    memcpy(data,&(e->data),bytes);
  return true;
}

void Taggable::setTagData(TagHandle* tag, void const* data)
{
  int i = findEntry(container,tag,size);
  size_t const bytes = tag->getBytes();
  Entry* e;
  if (i==-1)
  {
    e = addEntry(container,tag,size);
    if (bytes > sizeof(void*))
      e->data = malloc(bytes);
  }
  else
    e = container + i;
  if (bytes > sizeof(void*))
    memcpy(e->data,data,bytes);
  else
    memcpy(&(e->data),data,bytes);
}

const char* Taggable::getTagString(TagHandle* tag)
{
  assert(tag->getType()==PUMI_STR);
  int const i = findEntry(container,tag,size);
  if (i==-1)
    return 0;
  return static_cast<const char*>(container[i].data);
}

void Taggable::setTagString(TagHandle* tag, const char* data)
{
  assert(tag->getType()==PUMI_STR);
  int i = findEntry(container,tag,size);
  size_t const bytes = strlen(data)+1;
  Entry* e;
  if (i==-1)
  {
    e = addEntry(container,tag,size);
    e->data = 0;
  }
  else
    e = container + i;
  free(e->data);
  e->data = malloc(bytes);
  memcpy(e->data,data,bytes);
}

int Tag_GetType(pTag tag)
{
  return tag->getType();
}

int Tag_GetSize(pTag tag)
{
  return tag->getSize();
}

const char* Tag_GetName(pTag tag)
{
  return tag->getName().c_str();
}

/* getting around a fundamental design bug in std::set */
static pTag getTag(std::set<TagHandle>::iterator const& it)
{
  const TagHandle* p = &(*it);
  return const_cast<TagHandle*>(p);
}

pTag TagHolder_CreateTag(
    pTagHolder holder,
    const char* tag_name,
    int tag_type,
    int tag_size)
{
  TagHandle handle(tag_name,tag_type,tag_size);
  /* the handle will get copied, but just for this part its ok */
  return getTag(((holder->tags).insert(handle)).first);
}

int TagHolder_DelTag (pTagHolder holder, pTag tag)
{
  std::set<TagHandle>& tags = holder->tags;
  for (std::set<TagHandle>::iterator it = tags.begin();
      it != tags.end(); ++it)
    if (getTag(it) == tag)
    {
      tags.erase(it);
      break;
    }
  return PUMI_SUCCESS;
}

void TagHolder_ClearTag (pTagHolder holder)
{
  holder->tags.clear();
}

int TagHolder_CheckTag(pTagHolder holder, pTag tag, int tag_type)
{
  if (tag->getType() != tag_type)
    return PUMI_INVALID_TAG_HANDLE;
  std::set<TagHandle>& tags = holder->tags;
  for (std::set<TagHandle>::iterator it = tags.begin();
      it != tags.end(); ++it)
    if ( getTag(it) == tag )
      return PUMI_SUCCESS;
  return PUMI_INVALID_TAG_HANDLE;
}

int TagHolder_FindTag(pTagHolder holder, const char* tag_name, pTag& tag)
{
  TagHandle key(tag_name,0,0);
  std::set<TagHandle>::iterator it = holder->tags.find(key);
  if (it == holder->tags.end())
    return PUMI_INVALID_TAG_HANDLE;
  tag = getTag(it);
  return PUMI_SUCCESS;
}

int TagHolder_HasTag (pTagHolder holder, const pTag tag, int *exist)
{
  std::set<TagHandle>& tags = holder->tags;
  *exist = 0;
  for (std::set<TagHandle>::iterator it = tags.begin();
      it != tags.end(); ++it)
    if ( getTag(it) == tag )
      *exist = 1;
  return PUMI_SUCCESS;
}

void TagHolder_GetTag (pTagHolder holder, std::vector<pTag> &tags)
{
  std::set<TagHandle>& t = holder->tags;
  tags.resize(t.size());
  size_t i=0;
  for (std::set<TagHandle>::iterator it = t.begin();
      it != t.end(); ++it)
    tags[i++] = getTag(it);
}

bool Taggable_HasTag (pTaggable obj, pTag tag)
{
  /* this behavior is absolutely ridiculous, but we will
     keep it around for compatibility.... */
  if ( ! tag) return (obj->size != 0);
  return obj->hasTagData(tag);
}

void Taggable_DelTag (pTaggable obj, pTag tag)
{
  obj->deleteTagData(tag);
}

void Taggable_GetTag (pTaggable obj, std::vector<pTag>& tags)
{
  tags.resize(obj->size);
  for (int i=0; i < obj->size; ++i)
    tags[i] = obj->container[i].handle;
}
