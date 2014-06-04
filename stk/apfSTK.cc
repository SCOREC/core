/*
 * Copyright 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfSTK.h"
#include <apfMesh.h>
#include <apfShape.h>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

namespace apf {

/**
 *   \brief Implement an shards::ArrayDimTag for Quadrature points
 *
 *   Note that QPDimTag::Size does not dictate the size of that dimension,
 *   put_field does.
 */
struct QPDimTag : public shards::ArrayDimTag {
  enum { Size = 1 };                    ///< default size
  const char * name() const
  {
    return "apf::QPDimTag";
  }
  static const QPDimTag & tag()        ///< Singleton
  {
    static QPDimTag self;
    return self;
  }
};

typedef
stk::mesh::Field<double, stk::mesh::Cartesian, stk::mesh::Cartesian>
StkTensorField;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> StkVectorField;
typedef stk::mesh::Field<double> StkScalarField;
typedef
stk::mesh::Field<double, QPDimTag, stk::mesh::Cartesian, stk::mesh::Cartesian>
StkQPTensorField;
typedef
stk::mesh::Field<double, QPDimTag, stk::mesh::Cartesian >
StkQPVectorField;
typedef stk::mesh::Field<double, QPDimTag> StkQPScalarField;

typedef std::map<long,Node> GlobalMap;

template <class T>
T* makeStkField(
    const char* name,
    StkMetaData* metaData);

template <class T>
T* makeStkQPField(
    const char* name,
    const unsigned nqp,
    StkMetaData* metaData);

template <>
StkScalarField* makeStkField<StkScalarField>(
    const char* name,
    StkMetaData* metaData)
{
  StkScalarField* result;
  result = &(metaData->declare_field<StkScalarField>(name));
  stk::mesh::put_field(
      *result,
      metaData->node_rank(),
      metaData->universal_part());
  stk::io::set_field_role(*result,Ioss::Field::TRANSIENT);
  return result;
}

template <>
StkVectorField* makeStkField<StkVectorField>(
    const char* name,
    StkMetaData* metaData)
{
  StkVectorField* result;
  result = &(metaData->declare_field<StkVectorField>(name));
  stk::mesh::put_field(
      *result,
      metaData->node_rank(),
      metaData->universal_part(),
      3);
  stk::io::set_field_role(*result,Ioss::Field::TRANSIENT);
  return result;
}

template <>
StkTensorField* makeStkField<StkTensorField>(
    const char* name,
    StkMetaData* metaData)
{
  StkTensorField* result;
  result = &(metaData->declare_field<StkTensorField>(name));
  stk::mesh::put_field(
      *result,
      metaData->node_rank(),
      metaData->universal_part(),
      3,3);
  stk::io::set_field_role(*result,Ioss::Field::TRANSIENT);
  return result;
}

template <>
StkQPScalarField* makeStkQPField<StkQPScalarField>(
    const char* name,
    unsigned nqp,
    StkMetaData* metaData)
{
  StkQPScalarField* result;
  result = &(metaData->declare_field<StkQPScalarField>(name));
  stk::mesh::put_field(
      *result,
      metaData->element_rank(),
      metaData->universal_part(),
      nqp);
  stk::io::set_field_role(*result,Ioss::Field::TRANSIENT);
  return result;
}

template <>
StkQPVectorField* makeStkQPField<StkQPVectorField>(
    const char* name,
    unsigned nqp,
    StkMetaData* metaData)
{
  StkQPVectorField* result;
  result = &(metaData->declare_field<StkQPVectorField>(name));
  stk::mesh::put_field(
      *result,
      metaData->element_rank(),
      metaData->universal_part(),
      3,nqp);
  stk::io::set_field_role(*result,Ioss::Field::TRANSIENT);
  return result;
}

template <>
StkQPTensorField* makeStkQPField<StkQPTensorField>(
    const char* name,
    const unsigned nqp,
    StkMetaData* metaData)
{
  StkQPTensorField* result;
  result = &(metaData->declare_field<StkQPTensorField>(name));
  stk::mesh::put_field(
      *result,
      metaData->element_rank(),
      metaData->universal_part(),
      3,3,nqp);
  stk::io::set_field_role(*result,Ioss::Field::TRANSIENT);
  return result;
}

static Node lookup(long id, GlobalMap& map)
{
  assert(map.count(id));
  return map[id];
}

void writeStkField(
    Field* field,
    StkScalarField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToNodes)
{
  stk::mesh::BucketArray<StkScalarField> stkArray(*stkField,*bucket);
  int nodesInBucket = stkArray.dimension(0);
  for (int i=0; i < nodesInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    Node node = lookup(globalId,globalIdsToNodes);
    double value;
    value = getScalar(field, node.entity, node.node);
    stkArray(i) = value;
  }
}

void writeStkField(
    Field* field,
    StkVectorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToNodes)
{
  stk::mesh::BucketArray<StkVectorField> stkArray(*stkField,*bucket);
  int nodesInBucket = stkArray.dimension(1);
  for (int i=0; i < nodesInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    Node node = lookup(globalId,globalIdsToNodes);
    Vector3 value;
    getVector(field, node.entity, node.node, value);
    for (int j=0; j < 3; ++j)
      stkArray(j,i) = value[j];
  }
}

void writeStkField(
    Field* field,
    StkTensorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToNodes)
{
  stk::mesh::BucketArray<StkTensorField> stkArray(*stkField,*bucket);
  int nodesInBucket = stkArray.dimension(2);
  for (int i=0; i < nodesInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    Node node = lookup(globalId,globalIdsToNodes);
    Matrix3x3 value;
    getMatrix(field, node.entity, node.node, value);
    for (int j=0; j < 3; ++j)
    for (int k=0; k < 3; ++k)
      stkArray(k,j,i) = value[j][k];
  }
}

void writeStkField(
    Field* field,
    StkQPScalarField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToElems)
{
  stk::mesh::BucketArray<StkQPScalarField> stkArray(*stkField,*bucket);
  int nqp = stkArray.dimension(0);
  int elemsInBucket = stkArray.dimension(1);
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    for (int j=0; j < nqp; ++j)
      stkArray(j,i) = getScalar(field,elem,j);
  }
}

void writeStkField(
    Field* field,
    StkQPVectorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToElems)
{
  stk::mesh::BucketArray<StkQPVectorField> stkArray(*stkField,*bucket);
  int nqp = stkArray.dimension(1);
  int elemsInBucket = stkArray.dimension(2);
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Vector3 value;
    for (int j=0; j < nqp; ++j)
    {
      getVector(field,elem,j,value);
      for (int k=0; k < 3; ++k)
        stkArray(k,j,i) = value[k];
    }
  }
}

void writeStkField(
    Field* field,
    StkQPTensorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToElems)
{
  stk::mesh::BucketArray<StkQPTensorField> stkArray(*stkField,*bucket);
  int nqp = stkArray.dimension(2);
  int elemsInBucket = stkArray.dimension(3);
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Matrix3x3 value;
    for (int j=0; j < nqp; ++j)
    {
      getMatrix(field,elem,j,value);
      for (int k=0; k < 3; ++k)
      for (int l=0; l < 3; ++l)
        stkArray(l,k,j,i) = value[k][l];
    }
  }
}

void readStkField(
    Field* field,
    StkScalarField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToNodes)
{
  stk::mesh::BucketArray<StkScalarField> stkArray(*stkField,*bucket);
  int nodesInBucket = stkArray.dimension(0);
  for (int i=0; i < nodesInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    Node node = lookup(globalId,globalIdsToNodes);
    double value;
    value = stkArray(i);
    setScalar(field, node.entity, node.node, value);
  }
}

void readStkField(
    Field* field,
    StkVectorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToNodes)
{
  stk::mesh::BucketArray<StkVectorField> stkArray(*stkField,*bucket);
  int nodesInBucket = stkArray.dimension(1);
  for (int i=0; i < nodesInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    Node node = lookup(globalId,globalIdsToNodes);
    Vector3 value;
    for (int j=0; j < 3; ++j)
      value[j] = stkArray(j,i);
    setVector(field, node.entity, node.node, value);
  }
}

void readStkField(
    Field* field,
    StkTensorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToNodes)
{
  stk::mesh::BucketArray<StkTensorField> stkArray(*stkField,*bucket);
  int nodesInBucket = stkArray.dimension(2);
  for (int i=0; i < nodesInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    Node node = lookup(globalId,globalIdsToNodes);
    Matrix3x3 value;
    for (int j=0; j < 3; ++j)
    for (int k=0; k < 3; ++k)
      value[j][k] = stkArray(k,j,i);
    setMatrix(field, node.entity, node.node, value);
  }
}

void readStkField(
    Field* field,
    StkQPScalarField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToElems)
{
  stk::mesh::BucketArray<StkQPScalarField> stkArray(*stkField,*bucket);
  int nqp = stkArray.dimension(0);
  int elemsInBucket = stkArray.dimension(1);
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    for (int j=0; j < nqp; ++j)
      setScalar(field,elem,j,stkArray(j,i));
  }
}

void readStkField(
    Field* field,
    StkQPVectorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToElems)
{
  stk::mesh::BucketArray<StkQPVectorField> stkArray(*stkField,*bucket);
  int nqp = stkArray.dimension(1);
  int elemsInBucket = stkArray.dimension(2);
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Vector3 value;
    for (int j=0; j < nqp; ++j)
    {
      for (int k=0; k < 3; ++k)
        value[k] = stkArray(k,j,i);
      setVector(field,elem,j,value);
    }
  }
}

void readStkField(
    Field* field,
    StkQPTensorField* stkField,
    StkBucket* bucket,
    GlobalMap& globalIdsToElems)
{
  stk::mesh::BucketArray<StkQPTensorField> stkArray(*stkField,*bucket);
  int nqp = stkArray.dimension(2);
  int elemsInBucket = stkArray.dimension(3);
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = (*bucket)[i].identifier();
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    for (int j=0; j < nqp; ++j)
    {
      Matrix3x3 value;
      for (int k=0; k < 3; ++k)
      for (int l=0; l < 3; ++l)
        value[k][l] = stkArray(l,k,j,i);
      setMatrix(field,elem,j,value);
    }
  }
}

class StkBridge
{
  public:
    static StkBridge* get(
        Field* f,
        StkMetaData* metaData,
        bool exists);
    virtual ~StkBridge() {}
    void transfer(
        StkMetaData* metaData,
        StkBulkData* bulkData,
        GlobalMap& globalIdsToVerts,
        GlobalMap& globalIdsToElems,
        bool toStk)
    {
      stk::mesh::Selector overlapSelector =
        metaData->locally_owned_part() |
        metaData->globally_shared_part();
      stk::mesh::BucketVector buckets;
      stk::mesh::EntityRank rank = metaData->node_rank();
      if (isQP) rank = metaData->element_rank();
      stk::mesh::get_buckets(
          overlapSelector,
          bulkData->buckets(rank),
          buckets);
      GlobalMap* globalIdsToEnts = &globalIdsToVerts;
      if (isQP) globalIdsToEnts = &globalIdsToElems;
      APF_ITERATE(stk::mesh::BucketVector,buckets,it)
      {
        if (toStk)
          write(*it,*globalIdsToEnts);
        else
          read(*it,*globalIdsToEnts);
      }
    }
    virtual void read(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts) = 0;
    virtual void write(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts) = 0;
    Field* apfField;
    bool isQP;
};

template <class T>
class NodalBridge : public StkBridge
{
  public:
    NodalBridge(Field* f,
             StkMetaData* metaData,
             bool exists)
    {
      apfField = f;
      if (exists)
        stkField = metaData->get_field<T>(getName(f));
      else
        stkField = makeStkField<T>(getName(f),metaData);
      isQP = false;
    }
    virtual ~NodalBridge() {}
    virtual void read(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      readStkField(apfField,stkField,bucket,globalIdsToEnts);
    }
    virtual void write(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      writeStkField(apfField,stkField,bucket,globalIdsToEnts);
    }
  private:
    T* stkField;
};

template <class T>
class QPBridge : public StkBridge
{
  public:
    QPBridge(Field* f,
             StkMetaData* metaData,
             bool exists)
    {
      apfField = f;
      if (exists)
        stkField = metaData->get_field<T>(getName(f));
      else
      {
        FieldShape* shape = getShape(f);
        Mesh* m = getMesh(f);
        MeshIterator* it = m->begin(m->getDimension());
        MeshEntity* e = m->iterate(it);
        m->end(it);
        int nqp = shape->countNodesOn(m->getType(e));
        stkField = makeStkQPField<T>(getName(f),nqp,metaData);
      }
      isQP = true;
    }
    virtual ~QPBridge() {}
    virtual void read(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      readStkField(apfField,stkField,bucket,globalIdsToEnts);
    }
    virtual void write(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      writeStkField(apfField,stkField,bucket,globalIdsToEnts);
    }
  private:
    T* stkField;
};

StkBridge* StkBridge::get(
    Field* f,
    StkMetaData* metaData,
    bool exists)
{
  if (getShape(f) == getMesh(f)->getShape())
  {
    switch(getValueType(f))
    {
      case SCALAR:
        return new NodalBridge<StkScalarField>(f,metaData,exists);
      case VECTOR:
        return new NodalBridge<StkVectorField>(f,metaData,exists);
    //case MATRIX:
      default:
        return new NodalBridge<StkTensorField>(f,metaData,exists);
    }
  }
  else
  {
    switch(getValueType(f))
    {
      case SCALAR:
        return new QPBridge<StkQPScalarField>(f,metaData,exists);
      case VECTOR:
        return new QPBridge<StkQPVectorField>(f,metaData,exists);
    //case MATRIX:
      default:
        return new QPBridge<StkQPTensorField>(f,metaData,exists);
    }
  }
}

void generateGlobalIdsToEnts(
    GlobalNumbering* n,
    GlobalMap& globalIdsToEnts)
{
  DynamicArray<Node> nodes;
  getNodes(n, nodes);
  for (size_t i = 0; i < nodes.getSize(); ++i) {
    long id = getStkId(n, nodes[i]);
    globalIdsToEnts[id] = nodes[i];
  }
}

void declareField(Field* f, StkMetaData* md)
{
  delete StkBridge::get(f, md, false);
}

void copyFieldsToMeta(
    Mesh* m,
    StkMetaData* metaData)
{
  declareField(m->getCoordinateField(), metaData);
  for (int i=0; i < m->countFields(); ++i)
    declareField(m->getField(i), metaData);
}

void transferField(
    Field* f,
    StkMetaData* md,
    StkBulkData* bd,
    GlobalMap& vm,
    GlobalMap& em,
    bool toStk)
{
  StkBridge* bridge = StkBridge::get(f, md, true);
  bridge->transfer(md, bd, vm, em, toStk);
  delete bridge;
}

void transferFields(
    GlobalNumbering* n[4],
    StkMetaData* metaData,
    StkBulkData* bulkData,
    bool toStk)
{
  Mesh* m = getMesh(n[0]);
  GlobalMap nm;
  generateGlobalIdsToEnts(n[0], nm);
  GlobalMap em;
  generateGlobalIdsToEnts(n[m->getDimension()], em);
  transferField(m->getCoordinateField(), metaData, bulkData, nm, em, toStk);
  for (int i=0; i < m->countFields(); ++i)
    transferField(m->getField(i), metaData, bulkData, nm, em, toStk);
}

void copyFieldsToBulk(
    GlobalNumbering* n[4],
    StkMetaData* meta,
    StkBulkData* bulk)
{
  transferFields(n, meta, bulk, true);
}

void copyFieldsFromBulk(
    GlobalNumbering* n[4],
    StkMetaData* meta,
    StkBulkData* bulk)
{
  transferFields(n, meta, bulk, false);
}

const CellTopologyData* getTopology(Mesh* m, int t)
{
  FieldShape* s = m->getShape();
  if (t == Mesh::VERTEX)
    return shards::getCellTopologyData< shards::Node >();
/* right now this ignores  boundary layer entities, which we
   are unlikely to deal with in STK in the near future */
  if (s->getOrder()==1)
  {
    if (t == Mesh::EDGE)
      return shards::getCellTopologyData< shards::Line<2> >();
    if (t == Mesh::TRIANGLE)
      return shards::getCellTopologyData< shards::Triangle<3> >();
    if (t == Mesh::QUAD)
      return shards::getCellTopologyData< shards::Quadrilateral<4> >();
    if (t == Mesh::TET)
      return shards::getCellTopologyData< shards::Tetrahedron<4> >();
    if (t == Mesh::HEX)
      return shards::getCellTopologyData< shards::Hexahedron<8> >();
  }
/* we also assume that all second order bases being considered
   (lagrange and composite, so far) use these topologies */
  else if (s->getOrder()==2)
  {
    if (t == Mesh::EDGE)
      return shards::getCellTopologyData< shards::Line<3> >();
    if (t == Mesh::TRIANGLE)
      return shards::getCellTopologyData< shards::Triangle<6> >();
    if (t == Mesh::QUAD)
      return shards::getCellTopologyData< shards::Quadrilateral<8> >();
    if (t == Mesh::TET)
      return shards::getCellTopologyData< shards::Tetrahedron<10> >();
    if (t == Mesh::HEX)
      return shards::getCellTopologyData< shards::Hexahedron<20> >();
  }
  abort();
  return 0;
}

static int getFirstType(Mesh* m, int dim)
{
  MeshIterator* it = m->begin(dim);
  MeshEntity* e = m->iterate(it);
  m->end(it);
  return m->getType(e);
}

const CellTopologyData* getDimTopology(Mesh* m, int dim)
{
  return getTopology(m, getFirstType(m, dim));
}

const CellTopologyData* getCellTopology(Mesh* m)
{
  return getDimTopology(m, m->getDimension());
}

long getStkId(GlobalNumbering* numbers, Node node)
{
  return getNumber(numbers, node) + 1;
}

}
