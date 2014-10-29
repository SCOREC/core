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
#include <stk_mesh/base/FindRestriction.hpp>

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
  result = &(metaData->declare_field<StkScalarField>(
        stk::topology::NODE_RANK, name));
  stk::mesh::put_field(
      *result,
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
  result = &(metaData->declare_field<StkVectorField>(
        stk::topology::NODE_RANK, name));
  stk::mesh::put_field(
      *result,
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
  result = &(metaData->declare_field<StkTensorField>(
        stk::topology::NODE_RANK, name));
  stk::mesh::put_field(
      *result,
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
  result = &(metaData->declare_field<StkQPScalarField>(
        stk::topology::ELEMENT_RANK, name));
  stk::mesh::put_field(
      *result,
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
  result = &(metaData->declare_field<StkQPVectorField>(
        stk::topology::ELEMENT_RANK, name));
  stk::mesh::put_field(
      *result,
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
  result = &(metaData->declare_field<StkQPTensorField>(
        stk::topology::ELEMENT_RANK, name));
  stk::mesh::put_field(
      *result,
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
    StkScalarField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToNodes)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nodesInBucket = bucket.size();
  for (size_t i=0; i < nodesInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    Node node = lookup(globalId,globalIdsToNodes);
    double value;
    data[i] = getScalar(field, node.entity, node.node);
  }
}

void writeStkField(
    Field* field,
    StkVectorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToNodes)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nodesInBucket = bucket.size();
  for (size_t i=0; i < nodesInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    Node node = lookup(globalId,globalIdsToNodes);
    Vector3 value;
    getVector(field, node.entity, node.node, value);
    for (size_t j=0; j < 3; ++j)
      data[3*i + j] = value[j];
  }
}

void writeStkField(
    Field* field,
    StkTensorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToNodes)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nodesInBucket = bucket.size();
  for (size_t i=0; i < nodesInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    Node node = lookup(globalId,globalIdsToNodes);
    Matrix3x3 value;
    getMatrix(field, node.entity, node.node, value);
    for (size_t j=0; j < 3; ++j)
    for (size_t k=0; k < 3; ++k)
      data[9*i + 3*j + k] = value[j][k];
  }
}

void writeStkField(
    Field* field,
    StkQPScalarField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToElems)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nqp = stk::mesh::find_restriction(
      stkField, bucket.entity_rank(), bucket.supersets()).dimension();
  size_t elemsInBucket = bucket.size();
  for (size_t i=0; i < elemsInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    for (size_t j=0; j < nqp; ++j)
      data[nqp*i + j] = getScalar(field,elem,j);
  }
}

void writeStkField(
    Field* field,
    StkQPVectorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToElems)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nqp = stk::mesh::find_restriction(
      stkField, bucket.entity_rank(), bucket.supersets()).dimension();
  size_t elemsInBucket = bucket.size();
  for (size_t i=0; i < elemsInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Vector3 value;
    for (size_t j=0; j < nqp; ++j)
    {
      getVector(field,elem,j,value);
      for (size_t k=0; k < 3; ++k)
        data[nqp*3*i + 3*j + k] = value[k];
    }
  }
}

void writeStkField(
    Field* field,
    StkQPTensorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToElems)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nqp = stk::mesh::find_restriction(
      stkField, bucket.entity_rank(), bucket.supersets()).dimension();
  size_t elemsInBucket = bucket.size();
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Matrix3x3 value;
    for (int j=0; j < nqp; ++j)
    {
      getMatrix(field,elem,j,value);
      for (int k=0; k < 3; ++k)
      for (int l=0; l < 3; ++l)
        data[nqp*9*i + 9*j + 3*k + l] = value[k][l];
    }
  }
}

void readStkField(
    Field* field,
    StkScalarField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToNodes)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nodesInBucket = bucket.size();
  for (size_t i=0; i < nodesInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    Node node = lookup(globalId,globalIdsToNodes);
    double value;
    setScalar(field, node.entity, node.node, data[i]);
  }
}

void readStkField(
    Field* field,
    StkVectorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToNodes)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nodesInBucket = bucket.size();
  for (size_t i=0; i < nodesInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    Node node = lookup(globalId,globalIdsToNodes);
    Vector3 value;
    for (size_t j=0; j < 3; ++j)
      value[j] = data[3*i + j];
    setVector(field, node.entity, node.node, value);
  }
}

void readStkField(
    Field* field,
    StkTensorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToNodes)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nodesInBucket = bucket.size();
  for (size_t i=0; i < nodesInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    Node node = lookup(globalId,globalIdsToNodes);
    Matrix3x3 value;
    for (size_t j=0; j < 3; ++j)
    for (size_t k=0; k < 3; ++k)
      value[j][k] = data[9*i + 3*j + k];
    setMatrix(field, node.entity, node.node, value);
  }
}

void readStkField(
    Field* field,
    StkQPScalarField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToElems)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nqp = stk::mesh::find_restriction(
      stkField, bucket.entity_rank(), bucket.supersets()).dimension();
  size_t elemsInBucket = bucket.size();
  for (size_t i=0; i < elemsInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    for (size_t j=0; j < nqp; ++j)
      setScalar(field, elem, j, data[nqp*i + j]);
  }
}

void readStkField(
    Field* field,
    StkQPVectorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToElems)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nqp = stk::mesh::find_restriction(
      stkField, bucket.entity_rank(), bucket.supersets()).dimension();
  size_t elemsInBucket = bucket.size();
  for (size_t i=0; i < elemsInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Vector3 value;
    for (size_t j=0; j < nqp; ++j)
    {
      for (size_t k=0; k < 3; ++k)
        value[k] = data[nqp*3*i + 3*j + k];
      setVector(field,elem,j,value);
    }
  }
}

void readStkField(
    Field* field,
    StkQPTensorField& stkField,
    StkBucket& bucket,
    GlobalMap& globalIdsToElems)
{
  StkBulkData& bulk = bucket.mesh();
  double* data = stk::mesh::field_data(stkField, bucket);
  size_t nqp = stk::mesh::find_restriction(
      stkField, bucket.entity_rank(), bucket.supersets()).dimension();
  size_t elemsInBucket = bucket.size();
  for (int i=0; i < elemsInBucket; ++i)
  {
    long globalId = bulk.identifier(bucket[i]);
    MeshEntity* elem = lookup(globalId,globalIdsToElems).entity;
    Matrix3x3 value;
    for (int j=0; j < nqp; ++j)
    {
      getMatrix(field,elem,j,value);
      for (int k=0; k < 3; ++k)
      for (int l=0; l < 3; ++l)
        data[nqp*9*i + 9*j + 3*k + l] = value[k][l];
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
      stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
      if (isQP)
        rank = stk::topology::ELEMENT_RANK;
      bulkData->get_buckets(rank, overlapSelector, buckets);
      GlobalMap* globalIdsToEnts = &globalIdsToVerts;
      if (isQP)
        globalIdsToEnts = &globalIdsToElems;
      APF_CONST_ITERATE(stk::mesh::BucketVector, buckets, it)
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
        stkField = metaData->get_field<T>(stk::topology::NODE_RANK, getName(f));
      else
        stkField = makeStkField<T>(getName(f),metaData);
      isQP = false;
    }
    virtual ~NodalBridge() {}
    virtual void read(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      readStkField(apfField,*stkField,*bucket,globalIdsToEnts);
    }
    virtual void write(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      writeStkField(apfField,*stkField,*bucket,globalIdsToEnts);
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
        stkField = metaData->get_field<T>(stk::topology::ELEMENT_RANK, getName(f));
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
      readStkField(apfField,*stkField,*bucket,globalIdsToEnts);
    }
    virtual void write(
        StkBucket* bucket,
        GlobalMap& globalIdsToEnts)
    {
      writeStkField(apfField,*stkField,*bucket,globalIdsToEnts);
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
