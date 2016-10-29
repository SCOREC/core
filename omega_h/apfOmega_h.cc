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

static void field_from_osh(apf::Field* f, osh::Tag<osh::Real>* tag) {
  auto am = apf::getMesh(f);
  auto dim = am->getDimension();
  auto nc = tagbase->ncomps();
  auto data = osh::HostRead<osh::Real>(tag->array());
  auto vt = apf::getValueType(f);
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

static void field_from_osh(apf::Mesh* am, osh::Tag<osh::Real>* tag) {
  auto dim = am->getDimension();
  auto nc = tagbase->ncomps();
  int vt = -1;
  if (nc == 1) vt = apf::SCALAR;
  if (nc == dim) vt = apf::VECTOR;
  if (nc == dim * dim) vt = apf::MATRIX;
  auto f = apf::createGeneralField(am, tag->name().c_str(), vt, nc,
      am->getShape());
  field_from_osh(f, tag);
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

static void coords_to_osh(osh::Mesh* om, apf::Mesh* am) {
  field_to_osh(am->getCoordinateField(), om);
}

static void coords_from_osh(apf::Mesh* am, osh::Mesh* om) {
  field_from_osh(am->getCoordinateField(),
      om->get_tag<osh::Real>(0, "coordinates"));
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

static void
class_from_osh(apf::Mesh2* am, osh::Mesh2* om,
    std::vector<apf::MeshEntity*> const& ents,
    int ent_dim) {
{
  auto class_dim = osh::HostRead<osh::I8>(om->get_array<osh::I8>(ent_dim, "class_dim"));
  auto class_id = osh::HostRead<osh::LO>(om->get_array<osh::LO>(ent_dim, "class_id"));
  for (int i = 0; i < om->nents(ent_dim); ++i) {
    auto ge = am->findModelEntity(class_dim[i], class_id[i]);
    am->setModelEntity(ents[i], ge);
  }
}

static std::vector<apf::MeshEntity*>
verts_from_osh(osh::Mesh* om, apf::Mesh2* am) {
{
  std::vector<apf::MeshEntity*> verts(om->nverts());
  for (int i = 0; i < om->nverts(); ++i) {
    verts[i] = am->createVert(nullptr);
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

static void owners_from_osh(
    apf::Mesh2* am,
    osh::Mesh* om,
    std::vector<apf::MeshEntity*> const& ents,
    int ent_dim)
{
  auto owners = om->ask_owners(ent_dim);
  auto own_ranks = osh::HostRead<osh::I32>(owners.ranks);
  auto own_ids = osh::HostRead<osh::LO>(owners.idxs);
  if (om->parting() != OMEGA_H_ELEM_BASED) {
    /* currently MDS defines ownership between pairs of ranks
       only, which is insufficient to represent ghosting
       ownership of elements.
       therefore we will tag the owners to the entities,
       which can later be used by an apf::Sharing */
    apf::MeshTag* own_tag = am->findTag("owner");
    if (!own_tag) own_tag = am->createIntTag("owner", 1);
    for (int i = 0; i < om->nents(ent_dim); ++i) {
      am->setIntTag(ents[i], own_tag, &own_ranks[i]);
    }
  }
  PCU_Comm_Begin();
  for (int i = 0; i < osh_count(om, ent_dim); ++i) {
    PCU_COMM_PACK(own_ranks[i], own_ids[i]);
    PCU_COMM_PACK(own_ranks[i], ents[i]);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int own_id;
    PCU_COMM_UNPACK(own_id);
    apf::MeshEntity* r;
    PCU_COMM_UNPACK(r);
    int from = PCU_Comm_Sender();
    if (from == PCU_Comm_Self()) {
      assert((int)from == own_ranks[own_id]);
      continue;
    }
    apf::MeshEntity* owner = ents[own_id];
    am->addRemote(owner, from, r);
  }
  PCU_Comm_Begin();
  for (int i = 0; i < osh_count(om, ent_dim); ++i)
    if (own_ranks[i] == (int)(PCU_Comm_Self())) {
      apf::Copies remotes;
      am->getRemotes(ents[i], remotes);
      int ncopies = remotes.size();
      int self = PCU_Comm_Self();
      APF_ITERATE(apf::Copies, remotes, it) {
        PCU_COMM_PACK(it->first, it->second);
        PCU_COMM_PACK(it->first, ncopies);
        PCU_COMM_PACK(it->first, self);
        PCU_COMM_PACK(it->first, ents[i]);
        APF_ITERATE(apf::Copies, remotes, it2)
          if (it2->first != it->first) {
            PCU_COMM_PACK(it->first, it2->first);
            PCU_COMM_PACK(it->first, it2->second);
          }
      }
    }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    apf::MeshEntity* e;
    PCU_COMM_UNPACK(e);
    int ncopies;
    PCU_COMM_UNPACK(ncopies);
    for (int i = 0; i < ncopies; ++i) {
      int p;
      PCU_COMM_UNPACK(p);
      apf::MeshEntity* r;
      PCU_COMM_UNPACK(r);
      am->addRemote(e, p, r);
    }
  }
}

void from_omega_h(apf::Mesh2* am, osh::Mesh* om)
{
  std::vector<apf::MeshEntity*> ents[4];
  ents[0] = verts_from_osh(am, om);
  for (unsigned d = 1; d <= om->dim(); ++d)
    ents[d] = ents_from_osh(am, om, ents[0], d);
  coords_from_osh(am, om);
  for (unsigned d = 0; d <= om->dim(); ++d) {
    class_from_osh(am, om, ents[d], d);
    owners_from_osh(am, om, ents[d], d);
    apf::initResidence(am, d);
  }
  am->acceptChanges();
  fields_from_osh(am, om);
}

};
