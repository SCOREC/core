/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFNUMBERINGCLASS_H
#define APFNUMBERINGCLASS_H

#include "apfField.h"

namespace apf {

template <class T>
class NumberingOf : public FieldBase
{
  public:
    NumberingOf();
    virtual int countComponents() const;
    virtual int getScalarType();
    void init(const char* n,
              Mesh* m,
              FieldShape* s,
              int c);
    void init(Field* f);
    Field* getField();
    FieldDataOf<T>* getData();
    void getAll(MeshEntity* e, T* dat);
    T get(MeshEntity* e, int node, int component);
    void set(MeshEntity* e, int node, int component, T value);
  private:
    Field* field;
    int components;
};

}

#endif
