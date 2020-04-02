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

template <class T, class S = T>
class ElementOf : public Element
{
  public:
    ElementOf(FieldOf<S>* f, MeshEntity* e):
      Element(f,e)
    {
    }
    ElementOf(FieldOf<S>* f, VectorElement* p):
      Element(f,p)
    {
    }
    virtual ~ElementOf() {}
    S* getNodeValues()
    {
      return reinterpret_cast<S*>(&(this->nodeData[0]));
    }
    T getValue(Vector3 const& local)
    {
      T value[1];
      getComponents(local, reinterpret_cast<double*>(value));
      return value[0];
    }
    void getValues(NewArray<S>& values)
    {
      values.allocate(nen);
      S* nodeValues = getNodeValues();
      for (int i=0; i < nen; ++i)
        values[i] = nodeValues[i];
    }
};

}//namespace apf

#endif
