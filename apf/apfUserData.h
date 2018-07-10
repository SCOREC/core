/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_USER_DATA_H
#define APF_USER_DATA_H

#include "apfFieldData.h"

namespace apf
{
template <class T>
struct UserDataBase : public FieldDataOf<T>
{
  UserDataBase(FunctionBase<T> * f)
    : function(f)
  { }
  void init(FieldBase * f)
  {
    FieldDataOf<T>::field = f;
  }
  bool hasEntity(MeshEntity * e)
  {
    return FieldDataOf<T>::field->getShape()->countNodesOn(FieldDataOf<T>::field->getMesh()->getType(e)) > 0;
  }
  void removeEntity(MeshEntity * e) {}
  void get(MeshEntity * e, T * data)
  {
    function->eval(e,data);
  }
  void set(MeshEntity * e, T const * data) {}
  bool isFrozen() { return false; }
  virtual FieldData* clone()
  {
    FieldData* newData = new UserDataBase<T>(function);
    newData->init(FieldDataOf<T>::field);
    copyFieldData(static_cast<FieldDataOf<T>*>(newData),
                  static_cast<FieldDataOf<T>*>(FieldDataOf<T>::field->getData()));
    return newData;
  }
  FunctionBase<T> const* getFunction() const { return function; }
  void setFunction(FunctionBase<T>* func) { function = func; }
private:
  FunctionBase<T> * function;
};
typedef UserDataBase<double> UserData;
typedef UserDataBase<int> UserDataInt;
typedef UserDataBase<long> UserDataLong;
typedef UserDataBase<size_t> UserDataSizeT;
}

#endif
