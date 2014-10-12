/*
 * Copyright 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_STK_H
#define APF_STK_H

#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Shards_BasicTopologies.hpp>

namespace apf {

typedef stk::mesh::MetaData StkMetaData;
typedef stk::mesh::BulkData StkBulkData;
typedef stk::mesh::Bucket StkBucket;

struct StkModel
{
  int dim;
  int apfTag;
  std::string stkName;
};

typedef Array<DynamicArray<StkModel>, 4> StkModels;

void copyMeshToMeta(
    Mesh* m,
    StkModels& models,
    StkMetaData* meta);

void copyFieldsToMeta(
    Mesh* m,
    StkMetaData* meta);

void makeStkNumberings(Mesh* m, GlobalNumbering* n[4]);
void freeStkNumberings(Mesh* m, GlobalNumbering* n[4]);

void copyMeshToBulk(
    GlobalNumbering* n[4],
    StkModels& models,
    StkMetaData* meta,
    StkBulkData* bulk);

void copyFieldsToBulk(
    GlobalNumbering* n[4],
    StkMetaData* meta,
    StkBulkData* bulk);

void copyFieldsFromBulk(
    GlobalNumbering* n[4],
    StkMetaData* meta,
    StkBulkData* bulk);

const CellTopologyData* getDimTopology(Mesh* m, int dim);
const CellTopologyData* getCellTopology(Mesh* m);

int getLocalSideId(Mesh* m, MeshEntity* e,
    MeshEntity* side);

long getStkId(GlobalNumbering* numbers, Node node);

}

#endif
