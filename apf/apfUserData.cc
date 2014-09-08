/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfUserData.h"

namespace apf {

UserData::UserData(Function* f)
{
  function = f;
}

void UserData::init(FieldBase* f)
{
  field = f;
}

bool UserData::hasEntity(MeshEntity*)
{
  return true;
}

void UserData::removeEntity(MeshEntity*)
{
}

void UserData::get(MeshEntity* e, double* data)
{
  function->eval(e, data);
}

void UserData::set(MeshEntity*, double const*)
{
}

bool UserData::isFrozen()
{
  return false;
}

}
