/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFTAGDATA_H
#define APFTAGDATA_H

#include "apfFieldData.h"

namespace apf {

class TagMaker
{
  public:
    virtual MeshTag* make(Mesh* m, const char* n, int s) = 0;
};

template <class T>
class TagHelper;

template <>
class TagHelper<int> : public TagMaker
{
  public:
    MeshTag* make(Mesh* m, const char* n, int s)
    {
      return m->createIntTag(n,s);
    }
    void get(Mesh* m, MeshEntity* e, MeshTag* t, int* data)
    {
      return m->getIntTag(e,t,data);
    }
    void set(Mesh* m, MeshEntity* e, MeshTag* t, int const* data)
    {
      return m->setIntTag(e,t,data);
    }
};

template <>
class TagHelper<double> : public TagMaker
{
  public:
    MeshTag* make(Mesh* m, const char* n, int s)
    {
      return m->createDoubleTag(n,s);
    }
    void get(Mesh* m, MeshEntity* e, MeshTag* t, double* data)
    {
      return m->getDoubleTag(e,t,data);
    }
    void set(Mesh* m, MeshEntity* e, MeshTag* t, double const* data)
    {
      return m->setDoubleTag(e,t,data);
    }
};

template <>
class TagHelper<long> : public TagMaker
{
  public:
    MeshTag* make(Mesh* m, const char* n, int s)
    {
      return m->createLongTag(n,s);
    }
    void get(Mesh* m, MeshEntity* e, MeshTag* t, long* data)
    {
      return m->getLongTag(e,t,data);
    }
    void set(Mesh* m, MeshEntity* e, MeshTag* t, long const* data)
    {
      return m->setLongTag(e,t,data);
    }
};

class TagData
{
  public:
    void init(
        const char* name,
        Mesh* m,
        FieldShape* s,
        TagMaker* mk,
        int c);
    ~TagData();
    bool hasEntity(MeshEntity* e);
    void removeEntity(MeshEntity* e);
    MeshTag* getTag(MeshEntity* e);
    MeshTag* makeOrFindTag(const char* name, int size);
  private:
    void createTags(const char* name, int components);
    void detachTags();
    void destroyTags();
    Mesh* mesh;
    FieldShape* shape;
    TagMaker* maker;
    MeshTag* tags[Mesh::TYPES];
};

template <class T>
class TagDataOf : public FieldDataOf<T>
{
  public:
    virtual void init(FieldBase* f)
    {
      this->FieldData::field = f;
      mesh = f->getMesh();
      tagData.init(
                f->getName(),
                mesh,
                f->getShape(),
                &helper,
                f->countComponents());
    }
    virtual bool hasEntity(MeshEntity* e) {return tagData.hasEntity(e);}
    virtual void removeEntity(MeshEntity* e) {tagData.removeEntity(e);}
    virtual void get(MeshEntity* e, T* data)
    {
      helper.get(mesh,e,tagData.getTag(e),data);
    }
    virtual void set(MeshEntity* e, T const* data)
    {
      helper.set(mesh,e,tagData.getTag(e),data);
    }
    virtual bool isFrozen()
    {
      return false;
    }
    virtual FieldData* clone()
    {
      return new TagDataOf<T>();
    }
  private:
    Mesh* mesh;
    TagData tagData;
    TagHelper<T> helper;
};

}

#endif
