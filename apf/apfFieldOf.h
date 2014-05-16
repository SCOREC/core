/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFFIELDOF_H
#define APFFIELDOF_H

#include "apfField.h"
#include "apfFieldData.h"

namespace apf {

template <class T>
class FieldOf;

template <class T>
void project(FieldOf<T>* to, FieldOf<T>* from);

/* see Level 1 BLAS */
template <class T>
void axpy(double a, FieldOf<T>* x, FieldOf<T>* y);

template <class T>
class FieldOf : public Field
{
  public:
    virtual ~FieldOf() {}
    void setNodeValue(MeshEntity* e, int node, T const& value)
    {
      getData()->setNodeComponents(
          e,node,reinterpret_cast<double const*>(&value));
    }
    void getNodeValue(MeshEntity* e, int node, T& value)
    {
      getData()->getNodeComponents(
          e,node,reinterpret_cast<double*>(&value));
    }
    void project(Field* from)
    {
      apf::project<T>(this,static_cast<FieldOf<T>*>(from));
    }
    void axpy(double a, Field* x)
    {
      apf::axpy<T>(a,static_cast<FieldOf<T>*>(x),this);
    }
};

}

#endif
