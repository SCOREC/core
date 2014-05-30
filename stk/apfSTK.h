#ifndef APF_STK_H
#define APF_STK_H

#include <apf.h>
#include <apfMesh.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <Shards_BasicTopologies.hpp>

namespace apf {

typedef stk::mesh::fem::FEMMetaData StkMetaData;
typedef stk::mesh::BulkData StkBulkData;
typedef stk::mesh::Bucket StkBucket;

void copyToMetaData(
    Mesh* m,
    StkMetaData* metaData);

void copyToBulkData(
    Mesh* m,
    StkMetaData* metaData,
    StkBulkData* bulkData);

void copyFromBulkData(
    Mesh* m,
    StkMetaData* metaData,
    StkBulkData* bulkData);

const CellTopologyData* getCellTopology(Mesh* m);

}

#endif
