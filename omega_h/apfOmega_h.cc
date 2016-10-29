#include "apfOmega_h.h"

#include <vector>
#include <cassert>

#include <apfMesh2.h>
#include <apfNumbering.h>
#include <PCU.h>
#include <apf.h>
#include <Omega_h_math.hpp>

namespace apf {

static void components_to_osh(
    apf::Field* f,
    apf::MeshIterator* it,
    osh::HostWrite<osh::Real> data) {
  apf::MeshEntity* v;
  int i = 0;
  while ((v = am->iterate(it))) {
    apf::getComponents(f, v, 0, data.data() + i * nc);
    ++i;
  }
}

static void components_from_osh(
    apf::Field* f,
    apf::MeshIterator* it,
    osh::HostRead<osh::Real> data) {
  apf::MeshEntity* v;
  int i = 0;
  while ((v = am->iterate(it))) {
    apf::setComponents(f, v, 0, data.data() + i * nc);
    ++i;
  }
}

static void vectors_to_osh(
    apf::Field* f,
    apf::MeshIterator* it,
    osh::HostWrite<osh::Real> data) {
  auto dim = apf::getMesh(f)->getDimension();
  apf::MeshEntity* v;
  int i = 0;
  while ((v = am->iterate(it))) {
    apf::Vector3 x;
    apf::getVector(f, v, 0, x);
    for (int j = 0; j < dim; ++j) data[i * dim + j] = x[j];
    ++i;
  }
}

static void vectors_from_osh(
    apf::Field* f,
    apf::MeshIterator* it,
    osh::HostRead<osh::Real> data) {
  auto dim = apf::getMesh(f)->getDimension();
  apf::MeshEntity* v;
  int i = 0;
  while ((v = am->iterate(it))) {
    apf::Vector3 x(0,0,0);
    for (int j = 0; j < dim; ++j) x[j] = data[i * dim + j];
    apf::setVector(f, v, 0, x);
    ++i;
  }
}

template <int dim>
static void matrices_to_osh(
    apf::Field* f,
    apf::MeshIterator* it,
    osh::HostWrite<osh::Real> data) {
  apf::MeshEntity* v;
  int i = 0;
  while ((v = am->iterate(it))) {
    apf::Matrix3x3 x;
    apf::getMatrix(f, v, 0, x);
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
        data[k * dim + j] = x[j][k];
  }
}

template <int dim>
static void matrices_from_osh(
    apf::Field* f,
    apf::MeshIterator* it,
    osh::HostWrite<osh::Real> data) {
  apf::MeshEntity* v;
  int i = 0;
  while ((v = am->iterate(it))) {
    apf::Matrix3x3 x;
    for (int j = 0; j < dim; ++j)
      for (int k = 0; k < dim; ++k)
        x[j][k] = data[k * dim + j];
    apf::setMatrix(f, v, 0, x);
  }
}

static void field_to_osh(osh::Mesh* om, apf::Field* f) {
  auto dim = om->dim();
  auto am = apf::getMesh(f);
  std::string name = apf::getName(f);
  auto nc = apf::countComponents(f);
  auto vt = apf::getValueType(f);
  auto data = osh::HostWrite<osh::Real>(om->nverts() * nc);
  apf::MeshEntity* v;
  auto it = am->begin(0);
  int i = 0;
  if (vt == apf::VECTOR) {
    if (dim == 2) vectors_to_osh<2>(f, it, data);
    if (dim == 3) vectors_to_osh<3>(f, it, data);
  } if (vt == apf::MATRIX) {
    if (dim == 2) matrices_to_osh<2>(f, it, data);
    if (dim == 3) matrices_to_osh<3>(f, it, data);
  } else components_to_osh(f, it, data);
  auto xfer = OMEGA_H_LINEAR_INTERP;
  om->add_tag(osh::VERT, name, nc, xfer, OMEGA_H_DO_OUTPUT, data.write());
}

static void field_from_osh(apf::Mesh* am, osh::Tag<osh::Real>* tag) {
  auto dim = am->getDimension();
  auto nc = tagbase->ncomps();
  auto data = osh::HostRead<osh::Real>(tag->array());
  int vt = -1;
  if (nc == 1) vt = apf::SCALAR;
  if (nc == dim) vt = apf::VECTOR;
  if (nc == dim * dim) vt = apf::MATRIX;
  auto f = apf::createGeneralField(am, tag->name().c_str(), vt, nc,
      am->getShape());
  apf::MeshEntity* v;
  apf::MeshIterator* it = am->begin(0);
  if (vt == apf::VECTOR) {
    if (dim == 2) vectors_from_osh<2>(f, it, data);
    if (dim == 3) vectors_from_osh<3>(f, it, data);
  } if (vt == apf::MATRIX) {
    if (dim == 2) matrices_from_osh<2>(f, it, data);
    if (dim == 3) matrices_from_osh<3>(f, it, data);
  } else components_from_osh(f, it, data);
}

static void fields_to_osh(osh::Mesh* om, apf::Mesh* am) {
  for (int i = 0; i < am->countFields(); ++i)
    if (apf::getShape(am->getField(i)) == am->getShape())
      field_to_osh(om, am->getField(i));
}

static void fields_from_osh(apf::Mesh* am, osh::Mesh* om) {
  for (int i = 0; i < om->ntags(0); ++i)
    auto tagbase = om->get_tag(0, i);
    if (tagbase->type == OMEGA_H_F64 &&
        tagbase->name() != "coordinates" &&
        tagbase->name() != "param") {
      field_from_osh(am, dynamic_cast<Tag<osh::Real> const*>(tagbase));
    }
  }
}

static void coords_to_osh(apf::Mesh* mesh_apf, osh::Mesh* mesh_osh) {
  field_to_osh(mesh_apf->getCoordinateField(), mesh_osh);
}

static void class_to_osh(apf::Mesh* mesh_apf, osh::Mesh* mesh_osh, int dim) {
  auto nents = osh::LO(mesh_apf->count(dim));
  auto host_class_id = osh::HostWrite<osh::LO>(nents);
  auto host_class_dim = osh::HostWrite<osh::I8>(nents);
  auto iter = mesh_apf->begin(dim);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = mesh_apf->iterate(iter))) {
    auto me = mesh_apf->toModel(e);
    host_class_dim[i] = osh::I8(mesh_apf->getModelType(me));
    host_class_id[i] = mesh_apf->getModelTag(me);
    ++i;
  }
  mesh_apf->end(iter);
  mesh_osh->add_tag(dim, "class_dim", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT,
      osh::Read<osh::I8>(host_class_dim.write()));
  mesh_osh->add_tag(dim, "class_id", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT,
      osh::LOs(host_class_id.write()));
}

static void conn_to_osh(apf::Mesh* mesh_apf, osh::Mesh* mesh_osh,
    apf::Numbering* vert_nums, int d) {
  auto nhigh = osh::LO(mesh_apf->count(d));
  auto deg = d + 1;
  osh::HostWrite<osh::LO> host_ev2v(nhigh * deg);
  auto iter = mesh_apf->begin(d);
  apf::MeshEntity* he;
  int i = 0;
  while ((he = mesh_apf->iterate(iter))) {
    apf::Downward eev;
    auto deg2 = mesh_apf->getDownward(he, 0, eev);
    OMEGA_H_CHECK(deg == deg2);
    for (int j = 0; j < deg; ++j) {
      host_ev2v[i * deg + j] = apf::getNumber(vert_nums, eev[j], 0, 0);
    }
    ++i;
  }
  mesh_apf->end(iter);
  auto ev2v = osh::LOs(host_ev2v.write());
  osh::Adj high2low;
  if (d == 1) {
    high2low.ab2b = ev2v;
  } else {
    auto lv2v = mesh_osh->ask_verts_of(d - 1);
    auto v2l = mesh_osh->ask_up(0, d - 1);
    high2low = osh::reflect_down(ev2v, lv2v, v2l, d, d - 1);
  }
  mesh_osh->set_ents(d, high2low);
}

static void globals_to_osh(
    apf::Mesh* mesh_apf, osh::Mesh* mesh_osh, int dim) {
  apf::GlobalNumbering* globals_apf = apf::makeGlobal(
      apf::numberOwnedDimension(mesh_apf, "smb2osh_global", dim));
  apf::synchronize(globals_apf);
  auto nents = osh::LO(mesh_apf->count(dim));
  osh::HostWrite<osh::GO> host_globals(nents);
  auto iter = mesh_apf->begin(dim);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = mesh_apf->iterate(iter))) {
    host_globals[i++] = apf::getNumber(globals_apf, apf::Node(e, 0));
  }
  mesh_apf->end(iter);
  apf::destroyGlobalNumbering(globals_apf);
  auto globals = osh::Read<osh::GO>(host_globals.write());
  mesh_osh->add_tag(dim, "global", 1, OMEGA_H_GLOBAL, OMEGA_H_DO_OUTPUT,
      osh::Read<osh::GO>(host_globals.write()));
  auto owners = osh::owners_from_globals(
      mesh_osh->comm(), globals, osh::Read<osh::I32>());
  mesh_osh->set_owners(dim, owners);
}

void to_omega_h(apf::Mesh* mesh_apf, osh::Mesh* mesh_osh) {
  auto comm_mpi = PCU_Get_Comm();
  decltype(comm_mpi) comm_impl;
  MPI_Comm_dup(comm_mpi, &comm_impl);
  auto comm_osh = osh::CommPtr(new osh::Comm(comm_impl));
  mesh_osh->set_comm(comm_osh);
  auto dim = mesh_apf->getDimension();
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  mesh_osh->set_dim(mesh_apf->getDimension());
  mesh_osh->set_verts(osh::LO(mesh_apf->count(0)));
  coords_to_osh(mesh_apf, mesh_osh);
  class_to_osh(mesh_apf, mesh_osh, 0);
  globals_to_osh(mesh_apf, mesh_osh, 0);
  auto vert_nums = apf::numberOverlapDimension(mesh_apf, "apf2osh", 0);
  for (int d = 1; d <= dim; ++d) {
    conn_to_osh(mesh_apf, mesh_osh, vert_nums, d);
    class_to_osh(mesh_apf, mesh_osh, d);
    globals_to_osh(mesh_apf, mesh_osh, d);
  }
  apf::destroyNumbering(vert_nums);
  fields_to_osh(mesh_osh, mesh_apf);
}

static std::vector<apf::MeshEntity*>
verts_from_osh(osh::Mesh* om, apf::Mesh2* am) {
{
  std::vector<apf::MeshEntity*> verts(om->nverts());
  auto coords = osh::HostRead<osh::Real>(om->coords());
  osh::HostRead<osh::Real> param;
  if (om->has_tag(0, "param")) param = om->get_array<osh::Real>(0, "param");
  else param = osh::HostRead<osh::Real>(om->nverts() * 2, 0);
  auto class_dim = osh::HostRead<osh::I8>(om->get_array<osh::I8>(0, "class_dim"));
  auto class_id = osh::HostRead<osh::LO>(om->get_array<osh::LO>(0, "class_id"));
  for (int i = 0; i < om->nverts(); ++i) {
    auto ge = am->findModelEntity(class_dim[i], class_id[i]);
    apf::Vector3 x(coords.data() + i * 3);
    apf::Vector3 p(param.data() + i * 2);
    verts[i] = am->createVertex(ge, x, p);
  }
  assert(am->count(0) == om->nverts());
  return verts;
}

static std::vector<apf::MeshEntity*>
ents_from_osh(
    apf::Mesh2* am,
    osh::Mesh* om,
    std::vector<apf::MeshEntity*> const& verts,
    int ent_dim)
{
  std::vector<apf::MeshEntity*> ents(om->nents(ent_dim));
  auto ev2v = osh::HostRead<osh::LO>(om->ask_verts_of(ent_dim));
  auto class_dim = osh::HostRead<osh::I8>(om->get_array<osh::I8>(0, "class_dim"));
  auto class_id = osh::HostRead<osh::LO>(om->get_array<osh::LO>(0, "class_id"));
  apf::Mesh::Type t = apf::Mesh::simplexTypes[ent_dim];
  int nverts_per_ent = ent_dim + 1;
  for (int i = 0; i < osh_count(om, ent_dim); ++i) {
    auto ge = am->findModelEntity(class_dim[i], class_id[i]);
    apf::Downward ev;
    for (int j = 0; j < nverts_per_ent; ++j)
      ev[j] = verts[ev2v[i * nverts_per_ent + j]];
    ents[i] = apf::buildElement(am, ge, t, ev);
  }
  return ents;
}

};
