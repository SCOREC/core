/*
 * Copyright 2016 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <string>
#include <cassert>

#include "apfFile.h"
#include "apf.h"
#include "apfShape.h"
#include "apfTagData.h"

namespace apf {

static void save_string(FILE* file, const char* s) {
  fprintf(file, "%s\n", s);
}

static std::string restore_string(FILE* file) {
  std::string s;
  while (true) {
    int c = fgetc(file);
    assert(c != EOF);
    if (c == '\n') break;
    s.push_back(c);
  }
  return s;
}

static void save_int(FILE* file, int x) {
  fprintf(file, "%d\n", x);
}

static int restore_int(FILE* file) {
  int x;
  int ret = fscanf(file, "%d\n", &x);
  assert(ret == 1);
  return x;
}

static void save_field_meta(FILE* file, apf::Field* field) {
  save_string(file, getName(field));
  save_int(file, getValueType(field));
  save_int(file, countComponents(field));
  save_string(file, getShape(field)->getName());
}

static void restore_field_meta(FILE* file, apf::Mesh* mesh) {
  std::string field_name = restore_string(file);
  int value_type = restore_int(file);
  int ncomps = restore_int(file);
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  makeField(mesh, field_name.c_str(), value_type, ncomps,
      shape, new TagDataOf<double>);
}

void save_meta(FILE* file, apf::Mesh* mesh) {
  save_string(file, mesh->getShape()->getName());
  save_int(file, mesh->countFields());
  for (int i = 0; i < mesh->countFields(); ++i) {
    save_field_meta(file, mesh->getField(i));
  }
}

void restore_meta(FILE* file, apf::Mesh* mesh) {
  std::string shape_name = restore_string(file);
  apf::FieldShape* shape = getShapeByName(shape_name.c_str());
  if (shape != mesh->getShape()) mesh->changeShape(shape, false);
  int nfields = restore_int(file);
  for (int i = 0; i < nfields; ++i) {
    restore_field_meta(file, mesh);
  }
}

} // end namespace apf
