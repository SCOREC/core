#include "apfOmega_h.h"

#include <cassert>

#include <apfMesh2.h>
#include <apfNumbering.h>
#include <PCU.h>
#include <apf.h>

namespace apf {

static void coords_to_osh(apf::Mesh* mesh_apf, Omega_h::Mesh* mesh_osh) {
  auto dim = mesh_apf->getDimension();
  auto nverts = mesh_apf->count(0);
  Omega_h::HostWrite<Omega_h::Real> host_coords(Omega_h::LO(nverts) * dim);
  auto iter = mesh_apf->begin(0);
  apf::MeshEntity* v;
  int i = 0;
  while ((v = mesh_apf->iterate(iter))) {
    apf::Vector3 point;
    mesh_apf->getPoint(v, 0, point);
    for (int j = 0; j < dim; ++j) {
      host_coords[i * dim + j] = point[unsigned(j)];
    }
    ++i;
  }
  mesh_apf->end(iter);
  mesh_osh->add_tag(0, "coordinates", dim, OMEGA_H_LINEAR_INTERP,
      OMEGA_H_DO_OUTPUT, Omega_h::Reals(host_coords.write()));
}

static void class_to_osh(apf::Mesh* mesh_apf, Omega_h::Mesh* mesh_osh, int dim) {
  auto nents = Omega_h::LO(mesh_apf->count(dim));
  auto host_class_id = Omega_h::HostWrite<Omega_h::LO>(nents);
  auto host_class_dim = Omega_h::HostWrite<Omega_h::I8>(nents);
  auto iter = mesh_apf->begin(dim);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = mesh_apf->iterate(iter))) {
    auto me = mesh_apf->toModel(e);
    host_class_dim[i] = Omega_h::I8(mesh_apf->getModelType(me));
    host_class_id[i] = mesh_apf->getModelTag(me);
    ++i;
  }
  mesh_apf->end(iter);
  mesh_osh->add_tag(dim, "class_dim", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT,
      Omega_h::Read<Omega_h::I8>(host_class_dim.write()));
  mesh_osh->add_tag(dim, "class_id", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT,
      Omega_h::LOs(host_class_id.write()));
}

static void copy_conn(apf::Mesh* mesh_apf, Omega_h::Mesh* mesh_osh,
    apf::Numbering* vert_nums, int d) {
  auto nhigh = Omega_h::LO(mesh_apf->count(d));
  auto deg = d + 1;
  Omega_h::HostWrite<Omega_h::LO> host_ev2v(nhigh * deg);
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
  auto ev2v = Omega_h::LOs(host_ev2v.write());
  Omega_h::Adj high2low;
  if (d == 1) {
    high2low.ab2b = ev2v;
  } else {
    auto lv2v = mesh_osh->ask_verts_of(d - 1);
    auto v2l = mesh_osh->ask_up(0, d - 1);
    high2low = Omega_h::reflect_down(ev2v, lv2v, v2l, d, d - 1);
  }
  mesh_osh->set_ents(d, high2low);
}

static void globals_to_osh(
    apf::Mesh* mesh_apf, Omega_h::Mesh* mesh_osh, int dim) {
  apf::GlobalNumbering* globals_apf = apf::makeGlobal(
      apf::numberOwnedDimension(mesh_apf, "smb2osh_global", dim));
  apf::synchronize(globals_apf);
  auto nents = Omega_h::LO(mesh_apf->count(dim));
  Omega_h::HostWrite<Omega_h::GO> host_globals(nents);
  auto iter = mesh_apf->begin(dim);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = mesh_apf->iterate(iter))) {
    host_globals[i++] = apf::getNumber(globals_apf, apf::Node(e, 0));
  }
  mesh_apf->end(iter);
  apf::destroyGlobalNumbering(globals_apf);
  auto globals = Omega_h::Read<Omega_h::GO>(host_globals.write());
  mesh_osh->add_tag(dim, "global", 1, OMEGA_H_GLOBAL, OMEGA_H_DO_OUTPUT,
      Omega_h::Read<Omega_h::GO>(host_globals.write()));
  auto owners = Omega_h::owners_from_globals(
      mesh_osh->comm(), globals, Omega_h::Read<Omega_h::I32>());
  mesh_osh->set_owners(dim, owners);
}

static void apf2osh(apf::Mesh* mesh_apf, Omega_h::Mesh* mesh_osh) {
  auto comm_mpi = PCU_Get_Comm();
  decltype(comm_mpi) comm_impl;
  MPI_Comm_dup(comm_mpi, &comm_impl);
  auto comm_osh = Omega_h::CommPtr(new Omega_h::Comm(comm_impl));
  mesh_osh->set_comm(comm_osh);
  auto dim = mesh_apf->getDimension();
  CHECK(dim == 2 || dim == 3);
  mesh_osh->set_dim(mesh_apf->getDimension());
  mesh_osh->set_verts(Omega_h::LO(mesh_apf->count(0)));
  coords_to_osh(mesh_apf, mesh_osh);
  class_to_osh(mesh_apf, mesh_osh, 0);
  globals_to_osh(mesh_apf, mesh_osh, 0);
  auto vert_nums = apf::numberOverlapDimension(mesh_apf, "apf2osh", 0);
  for (int d = 1; d <= dim; ++d) {
    copy_conn(mesh_apf, mesh_osh, vert_nums, d);
    class_to_osh(mesh_apf, mesh_osh, d);
    globals_to_osh(mesh_apf, mesh_osh, d);
  }
  apf::destroyNumbering(vert_nums);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  if (argc != 4) {
    if (PCU_Comm_Self() == 0) {
      std::cout << "\n";
      std::cout << "usage: smb2osh in.dmg in.smb out.osh\n";
      std::cout << "   or: smb2osh                   (usage)\n";
    }
    PCU_Comm_Free();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_null();
  apf::Mesh2* am = apf::loadMdsMesh(argv[1], argv[2]);
  {
    auto lib = Omega_h::Library(&argc, &argv);
    auto world = lib.world();
    Omega_h::Mesh om(&lib);
    apf2osh(am, &om);
    am->destroyNative();
    apf::destroyMesh(am);
    Omega_h::binary::write(argv[3], &om);
  }
  PCU_Comm_Free();
  MPI_Finalize();
}

};
