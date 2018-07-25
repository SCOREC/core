/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFFIELDDATA_H
#define APFFIELDDATA_H

#include <string>
#include "apfField.h"
#include "apfShape.h"

namespace apf {

class FieldData
{
  public:
    virtual ~FieldData();
    virtual void init(FieldBase* f) = 0;
    virtual bool hasEntity(MeshEntity* e) = 0;
    virtual void removeEntity(MeshEntity* e) = 0;
    virtual bool isFrozen() = 0;
    virtual FieldData* clone();
    virtual void rename(const char* newName);
    FieldBase* getField() {return field;}
  protected:
    FieldBase* field;
};

template <class T>
class FieldDataOf;

// different reduction operations for Fields
template <class T>
class ReductionOp
{
  public:
    virtual T apply(T val1, T val2) const = 0;
};

template <class T>
class ReductionSum : public ReductionOp<T>
{
  T apply(T val1, T val2) const { return val1 + val2; };
};

template <class T>
class ReductionMin : public ReductionOp<T>
{
  T apply(T val1, T val2) const { return ( (val1 < val2) ? val1 : val2 ); };
};

template <class T>
class ReductionMax : public ReductionOp<T>
{
  T apply(T val1, T val2) const { return ( (val1 < val2) ? val2 : val1 ); };
};


/* instantiate (is this necessary with the global consts below?) */
template class ReductionSum<double>;
template class ReductionMin<double>;
template class ReductionMax<double>;



template <class T>
void synchronizeFieldData(FieldDataOf<T>* data, Sharing* shr, bool delete_shr=false);

void accumulateFieldData(FieldDataOf<double>* data, Sharing* shr, bool delete_shr=false);

void reduceFieldData(FieldDataOf<double>* data, Sharing* shr, bool delete_shr=false, const apf::ReductionOp<double>& reduce_op=ReductionSum<double>());

template <class T>
void copyFieldData(FieldDataOf<T>* from, FieldDataOf<T>* to);

template <class T>
void multiplyFieldData(FieldDataOf<T>* from, T d, FieldDataOf<T>* to);

template <class T>
void addFieldData(FieldDataOf<T>* from1, FieldDataOf<T>* from2, FieldDataOf<T>* to);

template <class T>
class FieldDataOf : public FieldData
{
  public:
    virtual void get(MeshEntity* e, T* data) = 0;
    virtual void set(MeshEntity* e, T const* data) = 0;
    void setNodeComponents(MeshEntity* e, int node, T const* components);
    void getNodeComponents(MeshEntity* e, int node, T* components);
    int getElementData(MeshEntity* entity, NewArray<T>& data);
};


} //namespace apf

#endif

