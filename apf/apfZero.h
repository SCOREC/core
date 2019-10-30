#ifndef APFZERO_H
#define APFZERO_H

#include "apfField.h"
#include "apfNew.h"

namespace apf
{

template <class T>
void setComponents(FieldBase* f, MeshEntity* e, int node, T const * components);

template <class T>
struct ZeroOp : public FieldOp
{
  ZeroOp(FieldBase * fb)
  {
    f = fb;
    int cmps = fb->countComponents();
    data.allocate(cmps);
    for (int i = 0; i < cmps; ++i)
      data[i] = 0;
    ent = 0;
  }
  bool inEntity(MeshEntity* e)
  {
    ent = e;
    return true;
  }
  void atNode(int node)
  {
    setComponents<T>(f,ent,node,&data[0]);
  }
  FieldBase* f;
  MeshEntity* ent;
  apf::NewArray<T> data;
};
}

#endif
