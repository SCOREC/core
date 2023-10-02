/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_USER_DATA_H
#define APF_USER_DATA_H

#include "apfFieldData.h"

namespace apf {

struct UserData : public FieldDataOf<double>
{
  UserData(Function* f);
  void init(FieldBase* f);
  bool hasEntity(MeshEntity* e);
  void removeEntity(MeshEntity* e);
  void get(MeshEntity* e, double* data);
  void set(MeshEntity* e, double const* data);
  bool isFrozen();
  virtual FieldData* clone();
  // using const * const gives an error on gcc/7.3.0 because the return is an
  // r-value which cannot be modified anyways
  Function const* getFunction() const { return function; }
  void setFunction(Function* func) { function = func; }
  private:
  Function* function;
};

}

#endif
