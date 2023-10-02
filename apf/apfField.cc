/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfField.h"
#include "apfShape.h"
#include "apfTagData.h"

namespace apf {

void FieldBase::init(
    const char* n,
    Mesh* m,
    FieldShape* s,
    FieldData* d)
{
  name = n;
  mesh = m;
  shape = s;
  data = d;
  d->init(this);
}

FieldBase::~FieldBase()
{
  delete data;
}

int FieldBase::countNodesOn(MeshEntity* e)
{
  return shape->countNodesOn(mesh->getType(e));
}

int FieldBase::countValuesOn(MeshEntity* e)
{
  return countNodesOn(e) * countComponents();
}

void FieldBase::changeData(FieldData* d)
{
  delete data;
  data = d;
}

void FieldBase::rename(const char* newName)
{
  data->rename(newName);
  name = newName;
}

FieldDataOf<double>* Field::getData()
{
  return static_cast<FieldDataOf<double>*>(data);
}

bool FieldOp::inEntity(MeshEntity*)
{
  return true;
}

void FieldOp::outEntity()
{
}

void FieldOp::atNode(int)
{
}

void FieldOp::apply(FieldBase* f)
{
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it)))
    {
      // this condition makes sure that the entity has node(s)
	  if ( s->countNodesOn(m->getType(e)) == 0)
		continue;
      if ( ! this->inEntity(e))
        continue;
      int n = f->countNodesOn(e);
      for (int i=0; i < n; ++i)
        this->atNode(i);
      this->outEntity();
    }
    m->end(it);
  }
}

struct ZeroOp : public FieldOp
{
  ZeroOp(Field* f)
  {
    field = f;
    int n = f->countComponents();
    data.allocate(n);
    for (int i = 0; i < n; ++i)
      data[i] = 0;
    ent = 0;
  }
  bool inEntity(MeshEntity* e)
  {
    ent = e;
    return true;
  }
  void atNode(int n)
  {
    setComponents(field, ent, n, &data[0]);
  }
  Field* field;
  MeshEntity* ent;
  apf::NewArray<double> data;
};

void zeroField(Field* f)
{
  ZeroOp op(f);
  op.apply(f);
}

} //namespace apf
