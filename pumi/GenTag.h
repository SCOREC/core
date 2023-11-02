/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef GENTAG_H
#define GENTAG_H

#include "pumi_errorcode.h"
#include <string>
#include <set>
#include <vector>

enum PUMI_TagType {/*0*/ PUMI_DBL, /*1*/ PUMI_INT, /*2*/ PUMI_LONG,
              /*3*/ PUMI_ENT,   /*4*/ PUMI_SET,  /*5*/ PUMI_PTR,
              /*6*/ PUMI_STR,   /*7*/ PUMI_BYTE, /*8*/ PUMI_TAGTYPES};

//******************************
class TagHandle
//******************************
{
  public:
    TagHandle(std::string const& in_name, const int in_type, const int in_size);
    ~TagHandle();
    std::string const& getName() {return tag_name;}
    int getType() {return tag_type;}
    int getSize() {return tag_size;}
    size_t getBytes() {return bytes;}
    bool operator<(const TagHandle& other) const;
    static size_t const typeSizes[PUMI_TAGTYPES];
  protected:
    std::string tag_name;
    int tag_type;
    int tag_size;
    size_t bytes;
};

//******************************
class TagHolder
//******************************
{
  public :
    std::set<TagHandle> tags; 
};

//******************************
class Taggable
//******************************
{
  public:
    struct Entry
    {
      TagHandle* handle;
      void* data;
    };
    Entry* container;
    int size;

    Taggable();
    ~Taggable();

    void deleteTagData(TagHandle* tag);
    void clearTagData();
    bool hasTagData (TagHandle* tag);
    bool getTagData(TagHandle* tag, void* data);
    void setTagData(TagHandle* tag, void const* data);
    const char* getTagString(TagHandle* tag);
    void setTagString(TagHandle* tag, const char* data);
};

// ***************************************************
//      Templated API for Tag
// ***************************************************

typedef class TagHandle* pTag;
typedef class Taggable* pTaggable;
typedef class TagHolder* pTagHolder;

int Tag_GetType (pTag tag);
int Tag_GetSize(pTag tag);
const char* Tag_GetName(pTag tag);

pTag TagHolder_CreateTag (pTagHolder holder, const char* tag_name, int tag_type, int tag_size);
int TagHolder_DelTag (pTagHolder holder, pTag tag);
void TagHolder_ClearTag (pTagHolder holder);
int TagHolder_CheckTag(pTagHolder holder, pTag tag, int tag_type);
int TagHolder_FindTag (pTagHolder holder, const char* tag_name, pTag& tag);
int TagHolder_HasTag (pTagHolder holder, const pTag tag, int *exist);
void TagHolder_GetTag (pTagHolder holder, std::vector<pTag> &tags);

bool Taggable_HasTag (pTaggable obj, pTag tag);
void Taggable_DelTag (pTaggable obj, pTag tag);
void Taggable_GetTag (pTaggable obj, std::vector<pTag>& tags);

template <typename T>
int Taggable_SetData (pTaggable obj, pTag tag, const T* data)
{
  obj->setTagData(tag,static_cast<void const*>(data));
  return PUMI_SUCCESS;
}

template <typename T>
int Taggable_GetData (pTaggable obj, pTag tag, T data[])
{
  if (obj->getTagData(tag,static_cast<void*>(data)))
    return PUMI_SUCCESS;
  else
    return PUMI_TAG_NOT_FOUND;
}

#endif  // _GENTAG_H_
