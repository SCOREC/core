/*
 * Copyright (C) 2011-2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "spr.h"
#include "apfMesh.h"

namespace spr {

apf::Field* getGradIPField(apf::Field* f, const char* name, int order)
{
  apf::Mesh* m = getMesh(f);
  int vt = apf::getValueType(f);
  assert(vt == apf::SCALAR || vt == apf::VECTOR);
  apf::Field* ip_field = apf::createIPField(m,name,vt+1,order);
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
  {
    apf::MeshElement* me = apf::createMeshElement(m,e);
    apf::Element* fe = apf::createElement(f,me);
    int np = countIntPoints(me,order);
    for (int p=0; p < np; ++p)
    {
      apf::Vector3 xi;
      apf::getIntPoint(me,order,p,xi);
      if (vt == apf::SCALAR)
      {
        apf::Vector3 value;
        apf::getGrad(fe,xi,value);
        apf::setVector(ip_field,e,p,value);
      }
      else
      {
        apf::Matrix3x3 value;
        apf::getVectorGrad(fe,xi,value);
        apf::setMatrix(ip_field,e,p,value);
      }
    }
    apf::destroyElement(fe);
    apf::destroyMeshElement(me);
  }
  m->end(it);
  return ip_field;
}

}
