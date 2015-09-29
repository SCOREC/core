/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apf_stkConfig.h>
#if !HAS_STK
#error "configuration bug"
#endif
#include "apfSTK.h"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ionit_Initializer.h>

namespace apf {

static void define_output_fields(stk::io::StkMeshIoBroker& mesh_data,
                          std::size_t output_file_idx)
{
  // Follow the approach in Albany::STKDiscretization::setupExodusOutput().
  const stk::mesh::MetaData& meta = mesh_data.meta_data();
  const stk::mesh::FieldVector& fields = meta.get_fields();
  for (std::size_t i = 0; i < fields.size(); i++) {
    try {
      mesh_data.add_field(output_file_idx, *fields[i]);
    }
    catch (std::runtime_error const&) {}
  }
}

void writeExodus(
    apf::Mesh* mesh,
    apf::StkModels& models,
    const char* filename,
    std::size_t& output_file_idx,
    const double time_val,
    Teuchos::RCP<stk::mesh::MetaData>& meta,
    Teuchos::RCP<stk::mesh::BulkData>& bulk,
    Teuchos::RCP<stk::io::StkMeshIoBroker>& mesh_data)
{
  apf::GlobalNumbering* n[4];
  apf::makeStkNumberings(mesh, n);

  if (meta.is_null()) {
    meta = Teuchos::rcp(new stk::mesh::MetaData(mesh->getDimension()));
    apf::copyMeshToMeta(mesh, models, meta.get());
    apf::copyFieldsToMeta(mesh, meta.get());
    meta->commit();
  }

  if (bulk.is_null()) {
    bulk = Teuchos::rcp(new stk::mesh::BulkData(*meta, MPI_COMM_WORLD));
    apf::copyMeshToBulk(n, models, meta.get(), bulk.get());
  }
  apf::copyFieldsToBulk(n, meta.get(), bulk.get());

  if (mesh_data.is_null()) {
    Ioss::Init::Initializer();
    mesh_data = Teuchos::rcp(new stk::io::StkMeshIoBroker(MPI_COMM_WORLD));
    output_file_idx = mesh_data->create_output_mesh(filename,
                                                    stk::io::WRITE_RESULTS);
    mesh_data->set_bulk_data(*bulk);
    define_output_fields(*mesh_data, output_file_idx);
  }

  mesh_data->process_output_request(output_file_idx, time_val);

  apf::freeStkNumberings(mesh, n);
}

}
