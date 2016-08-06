/*
 * Copyright 2016 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <string>
#include <cassert>
#include <cstdlib>

#include <reel.h>
#include <pcu_io.h>
#include "apfFile.h"
#include "apf.h"
#include "apfShape.h"
#include "apfTagData.h"

namespace apf {

static void save_string(pcu_file* file, const char* s) {
  pcu_write_string(file, s);
}

static std::string restore_string(pcu_file* file) {
  char* s;
  pcu_read_string(file, &s);
  std::string s_cpp(s);
  free(s);
  return s_cpp;
}

static void save_int(pcu_file* file, int x) {
  unsigned tmp = static_cast<unsigned>(x);
  PCU_WRITE_UNSIGNED(file, tmp);
}

static int restore_int(pcu_file* file) {
  unsigned tmp;
  PCU_READ_UNSIGNED(file, tmp);
  return static_cast<int>(tmp);
}

static void save_field_meta(pcu_file* file, apf::Field* field) {
  save_string(file, getName(field));
  save_int(file, getValueType(field));
  save_int(file, countComponents(field));
  save_string(file, getShape(field)->getName());
}

static void restore_field_meta(pcu_file* file, apf::Mesh* mesh) {
  std::string field_name = restore_string(file);
  int value_type = restore_int(file);
  int ncomps = restore_int(file);
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  if (shape == 0)
    reel_fail("field shape \"%s\" could not be found\n", shape_name.c_str());
  makeField(mesh, field_name.c_str(), value_type, ncomps,
      shape, new TagDataOf<double>);
}

static int latest_version_number = 2;

void save_meta(pcu_file* file, apf::Mesh* mesh) {
  save_string(file, mesh->getShape()->getName());
/* the first version of this system made the mistake of
 * not including a version number.
 * to introduce a backwards-compatible version number,
 * what we do here is to write a negative version
 * number where the old code would have expected the
 * number of fields.
 */
  save_int(file, -latest_version_number);
  save_int(file, mesh->countFields());
  for (int i = 0; i < mesh->countFields(); ++i) {
    save_field_meta(file, mesh->getField(i));
  }
}

void restore_meta(pcu_file* file, apf::Mesh* mesh) {
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  assert(shape != 0);
  if (shape != mesh->getShape()) mesh->changeShape(shape, false);
  int nfields_or_version = restore_int(file);
  int nfields;
  int version;
  if (nfields_or_version >= 0) {
    nfields = nfields_or_version;
    version = 1;
  } else {
    nfields = restore_int(file);
    version = -nfields_or_version;
  }
  assert(version <= latest_version_number);
  assert(nfields >= 0);
  assert(nfields < 256);
  for (int i = 0; i < nfields; ++i) {
    restore_field_meta(file, mesh);
  }
}

} // end namespace apf
