/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apf.h>
#include <dwrHierarchic.h>

namespace dwr {

apf::Field* createHierarchicField(apf::Mesh* m, const char* name,
    int valueType, int order)
{
  return createField(m,name,valueType,getHierarchic(order));
}

}
