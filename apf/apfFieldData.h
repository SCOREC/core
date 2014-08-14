/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFFIELDDATA_H
#define APFFIELDDATA_H

#include <string>
#include "apfField.h"
#include "apfShape.h"

namespace apf {

class FieldData
{
  public:
    virtual ~FieldData();
    virtual void init(FieldBase* f) = 0;
    virtual bool hasEntity(MeshEntity* e) = 0;
    virtual void removeEntity(MeshEntity* e) = 0;
    virtual void synchronize(Sharing* shr) = 0;
    virtual bool isFrozen() = 0;
    virtual FieldData* clone();
    FieldBase* getField() {return field;}
  protected:
    FieldBase* field;
};

template <class T>
class FieldDataOf;

template <class T>
void synchronizeFieldData(FieldDataOf<T>* data, Sharing* shr);

template <class T>
void copyFieldData(FieldDataOf<T>* to, FieldDataOf<T>* from);

void accumulateFieldData(FieldDataOf<double>* data, Sharing* shr);

template <class T>
class FieldDataOf : public FieldData
{
  public:
    virtual void synchronize(Sharing* shr)
    {
      synchronizeFieldData<T>(this, shr);
    }
    virtual void get(MeshEntity* e, T* data) = 0;
    virtual void set(MeshEntity* e, T const* data) = 0;
    void setNodeComponents(MeshEntity* e, int node, T const* components)
    {
      int n = field->countNodesOn(e);
      if (n==1)
        return set(e,components);
      int nc = field->countComponents();
      NewArray<T> allComponents(nc*n);
      if (this->hasEntity(e))
        get(e,&(allComponents[0]));
      for (int i=0; i < nc; ++i)
        allComponents[node*nc+i] = components[i];
      set(e,&(allComponents[0]));
    }
    void getNodeComponents(MeshEntity* e, int node, T* components)
    {
      int n = field->countNodesOn(e);
      if (n==1)
        return get(e,components);
      int nc = field->countComponents();
      NewArray<T> allComponents(nc*n);
      get(e,&(allComponents[0]));
      for (int i=0; i < nc; ++i)
        components[i] = allComponents[node*nc+i];
    }
    int getElementData(MeshEntity* entity, NewArray<T>& data)
    {
      Mesh* mesh = field->getMesh();
      int t = mesh->getType(entity);
      int ed = Mesh::typeDimension[t];
      FieldShape* fs = field->getShape();
      EntityShape* es = fs->getEntityShape(t);
      int nc = field->countComponents();
      int nen = es->countNodes();
      data.allocate(nc * nen);
      int n = 0;
      for (int d = 0; d <= ed; ++d)
      {
        if (fs->hasNodesIn(d))
        {
          Downward a;
          int na = mesh->getDownward(entity,d,a);
          for (int i = 0; i < na; ++i)
          {
            get(a[i],&(data[n]));
            n += nc * (fs->countNodesOn(mesh->getType(a[i])));
          }
        }
      }
      assert(n == nc * nen);
      return n;
    }
};

} //namespace apf

#endif

