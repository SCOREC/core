/*
 * Copyright 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_STK_H
#define APF_STK_H

#include "apfAlbany.h"
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Teuchos_RCP.hpp>

namespace apf {

typedef stk::mesh::MetaData StkMetaData;
typedef stk::mesh::BulkData StkBulkData;
typedef stk::mesh::Bucket StkBucket;

void copyMeshToMeta(
    Mesh* m,
    StkModels& models,
    StkMetaData* meta);

void copyFieldsToMeta(
    Mesh* m,
    StkMetaData* meta);

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

void writeExodus(
    apf::Mesh* mesh,
    apf::StkModels& models,
    const char* filename,
    const double time_val,
    Teuchos::RCP<stk::mesh::MetaData>& meta,
    Teuchos::RCP<stk::mesh::BulkData>& bulk,
    Teuchos::RCP<stk::io::StkMeshIoBroker>& mesh_data);

}

#endif
