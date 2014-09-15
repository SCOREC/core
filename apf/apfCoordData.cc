/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfCoordData.h"
#include "apf.h"
#include "apfMesh2.h"

namespace apf {

void CoordData::init(FieldBase* f)
{
  FieldData::field = f;
  mesh = f->getMesh();
  shape = f->getShape();
}

bool CoordData::hasEntity(MeshEntity* e)
{
  return FieldData::field->countNodesOn(e) != 0;
}

void CoordData::removeEntity(MeshEntity*)
{
  fail("apf::CoordData::removeEntity should not be called");
}

void CoordData::get(MeshEntity* e, double* data)
{
  Vector3* v = reinterpret_cast<Vector3*>(data);
  mesh->getPoint_(e,0,*v);
}

void CoordData::set(MeshEntity* e, double const* data)
{
  /* if the mesh is not a Mesh2 then this will break
     in an ugly way, but it seems to be the right
     place to have support for changing coordinates */
  Vector3 const* v = reinterpret_cast<Vector3 const*>(data);
  static_cast<Mesh2*>(mesh)->setPoint_(e,0,*v);
}

}

