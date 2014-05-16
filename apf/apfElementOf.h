/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFELEMENTOF_H
#define APFELEMENTOF_H

#include "apfElement.h"
#include "apfFieldOf.h"
#include "apfShape.h"

namespace apf {

template <class T>
class ElementOf : public Element
{
  public:
    ElementOf(FieldOf<T>* f, MeshEntity* e):
      Element(f,e)
    {
    }
    ElementOf(FieldOf<T>* f, VectorElement* p):
      Element(f,p)
    {
    }
    virtual ~ElementOf() {}
    T* getNodeValues()
    {
      return reinterpret_cast<T*>(&(this->nodeData[0]));
    }
    T getValue(Vector3 const& local)
    {
      T value;
      getComponents(local,reinterpret_cast<double*>(&value));
      return value;
    }
    void getValues(NewArray<T>& values)
    {
      values.allocate(nen);
      T* nodeValues = getNodeValues();
      for (int i=0; i < nen; ++i)
        values[i] = nodeValues[i];
    }
};

}//namespace apf

#endif
